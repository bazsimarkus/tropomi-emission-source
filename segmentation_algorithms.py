# -----------------------------------------------------------
# segmentation_algorithms.py
# Contains all algorithms and functions to open, handle, and classify TROPOMI SO2 products
# -----------------------------------------------------------

import math

import geopy.distance
import netCDF4
import numpy as np
from scipy import spatial
from skimage.morphology import flood_fill
from sklearn.cluster import DBSCAN

from CONFIG_PARAMETERS import molar_mass_SO2
from Volcano import *
from trajectory_algorithms import get_forward_trajectory, get_reverse_trajectory


def read_satellite_file(current_volcano):
    """Takes in a Volcano class, returns the arrays for SO2 density (PBL), lons, lats, and the SO2 detection mask"""

    print("Current volcano: " + current_volcano.name + " - " + str(current_volcano.time_begin) + 'T' + str(current_volcano.hour) + ':' + str(current_volcano.minute) + ':' + str(current_volcano.second) + 'Z')
    print(current_volcano.filename)

    # OPEN FILE
    fh = netCDF4.Dataset(current_volcano.filename, 'r', format="NETCDF4")

    global pixel_area_m2
    pixel_area_m2 = float(str(fh.__dict__['spatial_resolution'])[:-3].split('x')[0]) * float(str(fh.__dict__['spatial_resolution'])[:-3].split('x')[1]) * 1000000

    # GET DATA
    # - lat/lon
    lons = fh.groups['PRODUCT'].variables['longitude'][:][0, :, :]
    lats = fh.groups['PRODUCT'].variables['latitude'][:][0, :, :]

    # - SO2 maps
    SO2_PBL = fh.groups['PRODUCT'].variables['sulfurdioxide_total_vertical_column'][0, :, :]

    # - SO2 detection flag: 'sulfurdioxide_detection_flag'
    flag = fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables['sulfurdioxide_detection_flag'][0, :, :]
    mask = np.where(flag > 0, 1, 0)  # Assign 1 for flag values above 0

    SO2_1km = fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables['sulfurdioxide_total_vertical_column_1km'][0, :, :]
    SO2_7km = fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables['sulfurdioxide_total_vertical_column_7km'][0, :, :]
    SO2_15km = fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables['sulfurdioxide_total_vertical_column_15km'][0, :, :]

    return SO2_PBL, lons, lats, mask


def get_ground_truth(current_volcano, mask_shape):
    """Takes in a Volcano class, and the size of the detection mask (to fill up the rest with zeros), and returns the ground truth array"""

    # Ground truth from CSV
    gt_filename = ASSETS_DIRECTORY + "/ground_truth/" + str(os.path.basename(current_volcano.filename)) + ".csv"
    if os.path.exists(gt_filename):
        gt_data = np.genfromtxt(gt_filename, delimiter=',')
        gt_data = np.nan_to_num(gt_data)
        gt_shape = np.shape(gt_data)
        padded_array = np.zeros(mask_shape)
        padded_array[:gt_shape[0], :gt_shape[1]] = gt_data
        return padded_array
    else:
        return None


def get_lat_lon_index(latitude, longitude, lat_lon_combo):
    """Takes in decimal degrees latitude/longitude, and the combined lats/lons array, returns the index of the given arbitrary latitude/longitude in the combined array"""

    new_array = np.abs(lat_lon_combo - [latitude, longitude])
    flat_array = new_array.reshape(-1, 2)
    origin_coordinates = (0, 0)
    distance, index = spatial.KDTree(flat_array).query(origin_coordinates)
    latindex, lonindex = [int(index / 450), int(int(index) % 450)]

    return latindex, lonindex


def get_dbscan_label_map(SO2_PBL, mask):
    """Run the DBSCAN clustering algorithm on an array. Takes in the arrays for SO2 density (PBL), the SO2 detection mask, and returns an array of the detection mask, labelled with the generated cluster labels"""

    input_mask = np.array(mask.nonzero()).T  # nonzero - return the indeices of nonzero elements

    weight_matrix = []
    for current_point in input_mask:
        weight_matrix.append(SO2_PBL[current_point[0], current_point[1]]*1000)

    db = DBSCAN(eps=4, min_samples=3).fit(input_mask, sample_weight=weight_matrix)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    labels = labels + 1  # To eliminate the label "0", the minimum label should be 1, because we later filter with np.where

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)

    #print("Estimated number of clusters: %d" % n_clusters_)
    #print("Estimated number of noise points: %d" % n_noise_)

    mask_filled = np.where(mask > 0, 1, 0)  # Assign 1 for flag values above 0

    label_index = 0

    for current_index in input_mask:
        mask_filled[current_index[0]][current_index[1]] = labels[label_index]
        label_index = label_index + 1

    return mask_filled


def get_cluster_center_of_mass(SO2_PBL, mask_labelled, label):
    """Return the center of mass for a chosen cluster. Takes in the arrays for SO2 density (PBL), the DBSCAN-labeled SO2 detection mask, and the label of the chosen cluster, and returns the coordinates of the center of mass in the SO2 density array"""

    SO2_overzero = SO2_PBL
    SO2_overzero[SO2_overzero < 0] = 0
    SO2_overzero = np.power(SO2_overzero * 1000, 4)

    label_mask = np.where(mask_labelled == label, 1, 0)
    masses = np.where(label_mask > 0, SO2_overzero, 0)

    center_of_mass = np.around((masses * np.mgrid[0:masses.shape[0], 0:masses.shape[1]]).sum(1).sum(1) / masses.sum()).astype(int)

    return center_of_mass


def run_radius_100km(current_volcano): # Automatically associate every SO2 detection within a 100km radius (radius can be set)
    """Run the Radius search classifier (BC-1 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_circle = np.zeros(mask.shape)
    if mask.sum() > 0:
        radius_km = 100

        for current_lat_number in range(mask.shape[0]):  # Iterate through the mask and draw the circle with the specified radius
            for current_lon_number in range(mask.shape[1]):
                temp_diff = geopy.distance.geodesic((current_volcano.latitude, current_volcano.longitude), (lats[current_lat_number][current_lon_number], lons[current_lat_number][current_lon_number])).km
                if temp_diff < radius_km:
                    mask_circle[current_lat_number, current_lon_number] = 1

        mask_circle = np.multiply(mask_circle, mask) # Logical AND the circle with the original detection mask

    mask_filled = mask_circle
    mask_filled = np.where(mask_filled > 0, current_volcano.id, 0)
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 0, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "radius_100km"


def run_floodfill(current_volcano):
    """Run the Flood-fill classifier (BC-2 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_to_fill = mask
    if mask.sum() > 0:
        lat_lon_combo = np.dstack((lats, lons))
        latindex, lonindex = get_lat_lon_index(current_volcano.latitude, current_volcano.longitude, lat_lon_combo)

        if mask[latindex, lonindex] != 0:
            mask_to_fill = flood_fill(mask, (latindex, lonindex), 2, connectivity=2)
        else:  # Search for the nearest nonzero element
            search_area = 3
            idx = np.argwhere(mask)
            idx = idx[~(idx == [latindex, lonindex]).all(1)]
            nearest_index = idx[((idx - [latindex, lonindex]) ** 2).sum(1).argmin()]

            if abs(nearest_index[0] - latindex) < search_area and abs(nearest_index[1] - lonindex) < search_area:
                mask_to_fill = flood_fill(mask, (nearest_index[0], nearest_index[1]), 2, connectivity=2)

    mask_filled = np.where(mask_to_fill > 1, 1, 0)  # Assign 1 for flag values above 0
    mask_filled = np.where(mask_filled > 0, current_volcano.id, 0)
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 0, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "floodfill"


def run_floodfill_withwind(current_volcano):
    """Run the Radius search classifier with wind (BC-3 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_to_fill = mask
    trajectories = None
    if mask.sum() > 0:
        lat_lon_combo = np.dstack((lats, lons))
        latindex, lonindex = get_lat_lon_index(current_volcano.latitude, current_volcano.longitude, lat_lon_combo)

        if mask[latindex, lonindex] != 0:
            mask_to_fill = flood_fill(mask, (latindex, lonindex), 2, connectivity=2)
        else:  # Search for the nearest nonzero element
            search_area = 3
            idx = np.argwhere(mask)
            idx = idx[~(idx == [latindex, lonindex]).all(1)]
            nearest_index = idx[((idx - [latindex, lonindex]) ** 2).sum(1).argmin()]

            if abs(nearest_index[0] - latindex) < search_area and abs(nearest_index[1] - lonindex) < search_area:
                mask_to_fill = flood_fill(mask, (nearest_index[0], nearest_index[1]), 2, connectivity=2)

        trajectories = [get_forward_trajectory(current_volcano)]

        for index, current_latitude in enumerate(trajectories[0][1]):
            current_latindex, current_lonindex = get_lat_lon_index(current_latitude, trajectories[0][0][index], lat_lon_combo)
            if mask[current_latindex, current_lonindex] == 1:
                mask_to_fill = flood_fill(mask_to_fill, (current_latindex, current_lonindex), 2, connectivity=2)
            else:
                search_area = 3
                idx = np.argwhere(mask)
                idx = idx[~(idx == [current_latindex, current_lonindex]).all(1)]
                nearest_index = idx[((idx - [current_latindex, current_lonindex]) ** 2).sum(1).argmin()]

                if abs(nearest_index[0] - current_latindex) < search_area and abs(nearest_index[1] - current_lonindex) < search_area:
                    mask_to_fill = flood_fill(mask_to_fill, (nearest_index[0], nearest_index[1]), 2, connectivity=2)

    mask_filled = np.where(mask_to_fill > 1, 1, 0)  # Assign 1 for flag values above 0
    mask_filled = np.where(mask_filled > 0, current_volcano.id, 0)
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 0, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "floodfill_withwind", trajectories


def run_dbscan(current_volcano):
    """Run the DBSCAN classifier (BC-4 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_filled = np.zeros(mask.shape)
    if mask.sum() >= 4:     # Minimum number of samples in DBSCAN, under this the algorithm cannot run
        lat_lon_combo = np.dstack((lats, lons))
        latindex, lonindex = get_lat_lon_index(current_volcano.latitude, current_volcano.longitude, lat_lon_combo)

        mask_labelled = get_dbscan_label_map(SO2_PBL, mask)

        if mask_labelled[latindex, lonindex] != 0:
            value_to_fill = mask_labelled[latindex, lonindex]
        else:  # Search for the nearest nonzero element
            search_area = 3
            idx = np.argwhere(mask_labelled)
            idx = idx[~(idx == [latindex, lonindex]).all(1)]
            nearest_index = idx[((idx - [latindex, lonindex]) ** 2).sum(1).argmin()]

            if abs(nearest_index[0] - latindex) < search_area and abs(nearest_index[1] - lonindex) < search_area:
                value_to_fill = mask_labelled[nearest_index[0], nearest_index[1]]
            else:
                value_to_fill = -1

        mask_filled = np.where(mask_labelled == value_to_fill, 1, 0)  # Assign 1 for flag values above 0

    mask_filled = np.where(mask_filled > 0, current_volcano.id, 0)
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 0, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "dbscan"


def run_dbscan_withwind(current_volcano):
    """Run the DBSCAN classifier with wind (BC-5 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_filled = np.zeros(mask.shape)
    trajectories = None
    if mask.sum() >= 4:     # Minimum number of samples in DBSCAN, under this the algorithm cannot run
        lat_lon_combo = np.dstack((lats, lons))
        latindex, lonindex = get_lat_lon_index(current_volcano.latitude, current_volcano.longitude, lat_lon_combo)

        mask_labelled = get_dbscan_label_map(SO2_PBL, mask)

        if mask_labelled[latindex, lonindex] != 0:
            value_to_fill = mask_labelled[latindex, lonindex]
        else:  # Search for the nearest nonzero element
            search_area = 3
            idx = np.argwhere(mask_labelled)
            idx = idx[~(idx == [latindex, lonindex]).all(1)]
            nearest_index = idx[((idx - [latindex, lonindex]) ** 2).sum(1).argmin()]

            if abs(nearest_index[0] - latindex) < search_area and abs(nearest_index[1] - lonindex) < search_area:
                value_to_fill = mask_labelled[nearest_index[0], nearest_index[1]]
            else:
                value_to_fill = -1

        trajectories = [get_forward_trajectory(current_volcano)]

        for index, current_latitude in enumerate(trajectories[0][1]):
            current_latindex, current_lonindex = get_lat_lon_index(current_latitude, trajectories[0][0][index], lat_lon_combo)

            if mask[current_latindex, current_lonindex] == 1:
                mask_labelled = np.where(mask_labelled == mask_labelled[current_latindex, current_lonindex], value_to_fill, mask_labelled)
            else:
                search_area = 3
                idx = np.argwhere(mask)
                idx = idx[~(idx == [current_latindex, current_lonindex]).all(1)]
                nearest_index = idx[((idx - [current_latindex, current_lonindex]) ** 2).sum(1).argmin()]

                if abs(nearest_index[0] - current_latindex) < search_area and abs(nearest_index[1] - current_lonindex) < search_area:
                    mask_labelled = np.where(mask_labelled == mask_labelled[nearest_index[0], nearest_index[1]], value_to_fill,mask_labelled)

        mask_filled = np.where(mask_labelled == value_to_fill, 1, 0)  # Assign 1 for flag values above 0

    mask_filled = np.where(mask_filled > 0, current_volcano.id, 0)
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 0, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "dbscan_withwind", trajectories


def run_reverse_trajectory_dbscan(current_volcano): # Start reverse trajectory from center of mass point of cluster, if the trajectory is within range of the volcano, we associate the cluster
    """Run the Reverse Trajectory DBSCAN classifier (BC-6 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_to_fill = np.zeros(mask.shape)
    trajectories = None
    if mask.sum() >= 4:  # Minimum number of samples in DBSCAN, under this the algorithm cannot run
        mask_labelled = get_dbscan_label_map(SO2_PBL, mask)

        trajectories = []

        for unique_label in np.unique(mask_labelled):
            if unique_label != 0:
                cluster_center_of_mass = get_cluster_center_of_mass(SO2_PBL, mask_labelled, unique_label)
                height_of_cluster = current_volcano.elevation
                current_trajectory = get_reverse_trajectory(lats[cluster_center_of_mass[0], cluster_center_of_mass[1]], lons[cluster_center_of_mass[0], cluster_center_of_mass[1]], current_volcano.year, current_volcano.month, current_volcano.day, current_volcano.hour, height_of_cluster)
                trajectories.append(current_trajectory)

                for index, current_latitude in enumerate(current_trajectory[1]):
                    current_longitude = current_trajectory[0][index]
                    coordinate_distance = geopy.distance.geodesic((current_volcano.latitude, current_volcano.longitude),(current_latitude, current_longitude)).km

                    if coordinate_distance < 50: # If the trajectory is nearer than 50km, we associate the plume
                        label_mask = np.where(mask_labelled == unique_label, 1, 0)
                        mask_to_fill = np.logical_or(mask_to_fill, label_mask)
                        break   # If we find a volcano, we don't have to iterate further

    mask_filled = mask_to_fill
    mask_filled = np.where(mask_filled > 0, current_volcano.id, 0)
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 0, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "reverse_trajectory_dbscan", trajectories


def run_multi_dbscan(current_volcano):
    """Run the Multi-class DBSCAN classifier (MC-1 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)
    mask_to_fill = np.zeros(mask.shape)
    if mask.sum() >= 4:     # Minimum number of samples in DBSCAN, under this the algorithm cannot run
        volcano_centers = []
        for index, row in volcanolist_df.iterrows():
            volcano_centers.append([float(row['lat']), float(row['lon']), int(row['id'])])

        mask_labelled = get_dbscan_label_map(SO2_PBL, mask)
        mask_to_fill = mask_labelled

        # Build the list of clusters and their nearest clusters/volcanoes on the image
        clusters_on_image = []
        if len(np.unique(mask_labelled) > 1) and len(volcano_centers) > 0:
            for unique_label in np.unique(mask_labelled):
                cluster_response = {
                    'unique_label': None,
                    'center_of_mass_lon': None,
                    'center_of_mass_lat': None,
                    'nearest_cluster_label': None,
                    'nearest_cluster_distance': None,
                    'nearest_clusters': [],
                    'nearest_volcano_id': None,
                    'nearest_volcano_distance': None
                }

                nearest_volcano_id = None
                nearest_volcano_distance = None

                if unique_label != 0:
                    ccmx, ccmy = get_cluster_center_of_mass(SO2_PBL, mask_labelled, unique_label)
                    ccmlat, ccmlon = lats[ccmx, ccmy], lons[ccmx, ccmy]

                    for new_unique_label in np.unique(mask_labelled):
                        if new_unique_label != 0 and new_unique_label != unique_label:
                            nccmx, nccmy = get_cluster_center_of_mass(SO2_PBL, mask_labelled, new_unique_label)
                            nccmlat, nccmlon = lats[nccmx, nccmy], lons[nccmx, nccmy]

                            distance_to_new_cluster = geopy.distance.geodesic((ccmlat, ccmlon), (nccmlat, nccmlon)).km
                            cluster_response['nearest_clusters'].append([new_unique_label, distance_to_new_cluster])

                    cluster_response['nearest_clusters'].sort(key=lambda x: x[1])

                    nearest_cluster_label = None
                    nearest_cluster_distance = None

                    if len(cluster_response['nearest_clusters']) > 0:
                        nearest_cluster_label = cluster_response['nearest_clusters'][0][0]
                        nearest_cluster_distance = cluster_response['nearest_clusters'][0][1]

                    for current_volcano_center in volcano_centers:
                            if nearest_volcano_id is None:
                                nearest_volcano_id = current_volcano_center[2]
                                nearest_volcano_distance = geopy.distance.geodesic((ccmlat, ccmlon), (current_volcano_center[0], current_volcano_center[1])).km
                            else:
                                distance_to_new_volcano = geopy.distance.geodesic((ccmlat, ccmlon), (current_volcano_center[0], current_volcano_center[1])).km
                                if distance_to_new_volcano < nearest_volcano_distance:
                                    nearest_volcano_id = current_volcano_center[2]
                                    nearest_volcano_distance = distance_to_new_volcano

                    cluster_response['unique_label'] = unique_label
                    cluster_response['center_of_mass_lat'] = ccmlat
                    cluster_response['center_of_mass_lon'] = ccmlon
                    cluster_response['nearest_cluster_label'] = nearest_cluster_label
                    cluster_response['nearest_cluster_distance'] = nearest_cluster_distance
                    cluster_response['nearest_volcano_id'] = nearest_volcano_id
                    cluster_response['nearest_volcano_distance'] = nearest_volcano_distance

                    clusters_on_image.append(cluster_response)

            # Beginning of the algorithm

            tolerance_km = 200  # Tolerance parameter for the algorithm

            associations = []
            closest_cluster_to_volcano = sorted(clusters_on_image, key=lambda x: x['nearest_volcano_distance'])[0]

            previous_cluster = None
            current_cluster = closest_cluster_to_volcano
            source_volcano_id = closest_cluster_to_volcano['nearest_volcano_id'] # The original source volcano ID, this value gets replaced, when it's not the nearest volcano anymore

            shortest_index = 1 # We keep track in this variable, how many times we have switched source volcanoes, aka if we are jumping to the second, third, fourth... shortest volcano-cluster distance
            while True: # We start growing the associations list, and stop if we have associated every cluster
                if current_cluster['nearest_volcano_id'] == source_volcano_id:
                    if current_cluster['unique_label'] not in [i[0] for i in associations]:
                        associations.append([current_cluster['unique_label'], current_cluster['nearest_volcano_id']])
                    else:
                        for current_near_cluster in current_cluster['nearest_clusters']: # As the nearest_clusters list is already sorted, we can just iterate through it and stop when we find the first unassociated cluster
                            if current_near_cluster[0] not in [i[0] for i in associations]:
                                previous_cluster = current_cluster
                                current_cluster = list(filter(lambda t: t['unique_label'] == current_near_cluster[0], clusters_on_image))[0]
                                break
                else:
                    if current_cluster['nearest_volcano_distance'] > tolerance_km and current_cluster['nearest_volcano_distance'] > geopy.distance.geodesic((current_cluster['center_of_mass_lat'], current_cluster['center_of_mass_lon']), (previous_cluster['center_of_mass_lat'], previous_cluster['center_of_mass_lon'])).km:
                        associations.append([current_cluster['unique_label'], source_volcano_id])
                    else:
                        shortest_index = shortest_index + 1  # Increase the shortest_index for the future changes
                        closest_cluster_to_volcano = sorted(clusters_on_image, key=lambda x: x['nearest_volcano_distance'])[shortest_index - 1]  # Get the cluster with the n-th nearest volcano-cluster distance

                        previous_cluster = current_cluster
                        current_cluster = closest_cluster_to_volcano  # Change current cluster to the new cluster
                        source_volcano_id = closest_cluster_to_volcano['nearest_volcano_id']  # Change the source volcano to the new nearest volcano

                if len(associations) == (len(np.unique(mask_labelled)) - 1): # If we have associated every cluster, we stop the while loop
                    break

            # Remove clusters that are far away from any volcano or cluster from the association list, using the tolerance
            for current_associated_cluster in clusters_on_image:
                if current_associated_cluster['nearest_volcano_distance'] is not None and current_associated_cluster['nearest_cluster_distance'] is not None:
                    if current_associated_cluster['nearest_volcano_distance'] > tolerance_km and current_associated_cluster['nearest_cluster_distance'] > tolerance_km * 2:
                        for associations_to_remove in associations:
                            if associations_to_remove[0] == current_associated_cluster['unique_label']:
                                associations_to_remove[1] = -1

            for current_association in associations:
                mask_to_fill[mask_to_fill == current_association[0]] = current_association[1]

            #print(associations)

    mask_filled = mask_to_fill
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 1, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "multi_dbscan"


def run_multi_reverse_trajectory_dbscan(current_volcano): # Start reverse trajectory from center of mass point of cluster, if the trajectory is within range of the volcano, we associate the cluster
    """Run the Multi-class Reverse Trajectory DBSCAN classifier (MC-2 in Markus et al. 2022 Frontiers). Takes in a Volcano class, and returns all relevant information, including the filled detection mask with the volcano ID-s"""

    SO2_PBL, lons, lats, mask = read_satellite_file(current_volcano)
    gt_data = get_ground_truth(current_volcano, mask.shape)

    mask_to_fill = np.zeros(mask.shape)
    trajectories = None
    if mask.sum() >= 4:  # Minimum number of samples in DBSCAN, under this the algorithm cannot run
        image_center_point = [current_volcano.latitude, current_volcano.longitude]
        volcano_centers = []

        for index, row in volcanolist_df.iterrows():
            if abs(image_center_point[0] - float(row['lat'])) < 20 and abs(image_center_point[1] - float(row['lon'])) < 20: # Only append volanoes that are maximum +/- 20 degrees far from the center point
                volcano_centers.append([float(row['lat']), float(row['lon']), int(row['id'])])

        mask_labelled = get_dbscan_label_map(SO2_PBL, mask)
        mask_to_fill = mask_labelled

        associations = []

        if len(np.unique(mask_labelled) > 1) and len(volcano_centers) > 0:
            trajectories = []
            for unique_label in np.unique(mask_labelled):
                if unique_label != 0:
                    cluster_center_of_mass = get_cluster_center_of_mass(SO2_PBL, mask_labelled, unique_label)

                    nearest_volcano_id = None
                    nearest_volcano_distance = None

                    for current_volcano_center in volcano_centers:
                            if nearest_volcano_id is None:
                                nearest_volcano_id = current_volcano_center[2]
                                nearest_volcano_distance = geopy.distance.geodesic((lats[cluster_center_of_mass[0], cluster_center_of_mass[1]], lons[cluster_center_of_mass[0], cluster_center_of_mass[1]]), (current_volcano_center[0], current_volcano_center[1])).km
                            else:
                                distance_to_new_volcano = geopy.distance.geodesic((lats[cluster_center_of_mass[0], cluster_center_of_mass[1]], lons[cluster_center_of_mass[0], cluster_center_of_mass[1]]), (current_volcano_center[0], current_volcano_center[1])).km
                                if distance_to_new_volcano < nearest_volcano_distance:
                                    nearest_volcano_id = current_volcano_center[2]
                                    nearest_volcano_distance = distance_to_new_volcano

                    height_of_cluster = int(volcanolist_df.loc[volcanolist_df['id'] == nearest_volcano_id]['elevation'])

                    current_trajectory = get_reverse_trajectory(lats[cluster_center_of_mass[0], cluster_center_of_mass[1]], lons[cluster_center_of_mass[0], cluster_center_of_mass[1]], current_volcano.year, current_volcano.month, current_volcano.day, current_volcano.hour, height_of_cluster)
                    trajectories.append(current_trajectory)

                    nearest_volcano_id = None
                    nearest_volcano_distance = None

                    for index, current_latitude in enumerate(current_trajectory[1]):
                        current_longitude = current_trajectory[0][index]
                        for current_volcano_center in volcano_centers:
                            if nearest_volcano_id is None:
                                nearest_volcano_id = current_volcano_center[2]
                                nearest_volcano_distance = geopy.distance.geodesic((current_latitude, current_longitude), (current_volcano_center[0], current_volcano_center[1])).km
                            else:
                                distance_to_new_volcano = geopy.distance.geodesic((current_latitude, current_longitude), (current_volcano_center[0], current_volcano_center[1])).km
                                if distance_to_new_volcano < nearest_volcano_distance:
                                    nearest_volcano_id = current_volcano_center[2]
                                    nearest_volcano_distance = distance_to_new_volcano

                    if nearest_volcano_distance is not None and nearest_volcano_distance < 200:
                        associations.append([unique_label, nearest_volcano_id])
                    else:
                        associations.append([unique_label, -1])

            for current_association in associations:
                mask_to_fill[mask_to_fill == current_association[0]] = current_association[1]

        #print(associations)

    mask_filled = mask_to_fill
    center_volcano_mask = np.where(mask_filled == current_volcano.id, 1, 0)
    mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
    mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2

    return 1, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, "multi_reverse_trajectory_dbscan", trajectories
