# -----------------------------------------------------------
# plotting_functions.py
# Contains all functions to plot the results of function in segmentation_algorithms.py
# -----------------------------------------------------------

import matplotlib.colors
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.ticker import FixedLocator
from mpl_toolkits.basemap import Basemap

import segmentation_algorithms
from Volcano import *


def plot_segmentation_result(sem_inst, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, method, trajectories=None):
    """Plot the result of an algorithm. Takes in the result of an algorithm function in segmentation_algorithms.py, returns the plot figure"""

    PRINT_MODE = 0 # set to 1 if plotting for document, set to 0 if viewed on PC
    zoom_span = 6 # Zooming parameter on plot (the larger, the more zoomed out)

    # PLOT
    # Set plot DPI high for Spyder plot
    if PRINT_MODE == 1:
        plt.rcParams['figure.dpi'] = 500
        fig = plt.figure(figsize=(4, 2))
    else:
        plt.rcParams['figure.dpi'] = 200
        fig = plt.figure(figsize=(16, 9))

    m = Basemap(llcrnrlon=current_volcano.longitude - zoom_span / 2,
                llcrnrlat=current_volcano.latitude - zoom_span / 2,
                urcrnrlon=current_volcano.longitude + zoom_span / 2,
                urcrnrlat=current_volcano.latitude + zoom_span / 2, resolution='l', projection='cyl',
                lat_ts=40, lat_0=current_volcano.latitude, lon_0=current_volcano.longitude)

    m.drawparallels(np.arange(-85., 85., 2.), labels=[1, 0, 0, 0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 2.), labels=[0, 0, 0, 1], fontsize=10)
    try:
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
    except:
        print("No basemap lines drawn")

    # - SO2 map
    vmin, vmax = 0, 0.01
    vcd2plot, vcd2plot_str = SO2_PBL, 'SO2_PBL'
    cs = m.pcolormesh(lons, lats, vcd2plot, latlon=True, vmin=vmin, vmax=vmax, cmap='CMRmap_r', shading='auto')

    # - SO2 detection mask
    colorbar_labels = []

    if sem_inst == 0:
        fillcmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["w", "g"], 2)
        cmask_filled = m.pcolormesh(lons, lats, mask_filled.astype(np.float32), latlon=True, vmin=0, vmax=1, cmap=fillcmap, alpha=.5, shading='auto')
    else:
        # Crop the clusters that are not on the plot, to eliminate non-present colorbar labels
        for iterator_row in range(len(mask_filled)):
            for iterator_column in range(len(mask_filled[0])):
                if lons[iterator_row, iterator_column] > current_volcano.longitude + zoom_span / 2 or lons[iterator_row, iterator_column] < current_volcano.longitude - zoom_span / 2:
                    mask_filled[iterator_row, iterator_column] = 0
                if lats[iterator_row, iterator_column] > current_volcano.latitude + zoom_span / 2 or lats[iterator_row, iterator_column] < current_volcano.latitude - zoom_span / 2:
                    mask_filled[iterator_row, iterator_column] = 0

        # Count the remaining unique elements after cropping the image
        mask_filled = mask_filled.astype(int)
        number_of_remained_unique_elements_after_crop = len(np.unique(mask_filled))
        mask_filled = mask_filled.astype(np.float32)

        for element in np.unique(mask_filled[~np.isnan(mask_filled)]):
            current_id = int(element)
            if current_id != 0:
                current_label = int(current_id)
                colorbar_labels.append(current_label)
        if len(colorbar_labels) == 0:
            colorbar_labels.append(0)
            mycmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["w", "w"])
            cmask_filled = m.pcolormesh(lons, lats, mask_filled, latlon=True, cmap=mycmap, alpha=.8, shading='auto')
        else:
            mask_filled[mask_filled == 0] = np.nan  # to make the zeros transparent in pcolormesh
            mycmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["C0", "yellow", "C1", "C3", "C2"], number_of_remained_unique_elements_after_crop - 1)
            cmask_filled = m.pcolormesh(lons, lats, mask_filled, latlon=True, cmap=mycmap, alpha=.8, shading='auto')

    # - SO2 detection mask
    cmask = m.pcolormesh(lons, lats, mask.astype(np.float32), latlon=True, vmin=0, vmax=1, cmap=plt.cm.get_cmap('Greys', 2), alpha=.25, shading='auto')

    other_volcanoes = []
    # - cosmetics
    for index, row in volcanolist_df.iterrows():
        plt.text(row['lon'], row['lat'], str(row['fullname']) + " (" + str(row['id']) + ")", clip_on=True, fontsize=6)
        plt.plot(row['lon'], row['lat'], '^b', zorder=2)

        if (abs(float(row['lon']) - float(current_volcano.longitude)) < zoom_span / 2) and abs(
                (float(row['lat']) - float(current_volcano.latitude)) < zoom_span / 2) and str(
                row['name']) != current_volcano.name:
            other_volcanoes.append(row)

    if trajectories is not None:
        for current_trajectory in trajectories:
            plt.plot(current_trajectory[0], current_trajectory[1], color='C9', label='Plume trajectory')

        if PRINT_MODE == 1:
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), prop={'size': 4}, loc='upper right')
        else:
            handles, labels = plt.gca().get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            plt.legend(by_label.values(), by_label.keys(), prop={'size': 10}, loc='upper right')

    if PRINT_MODE == 0:
        title = '$SO_{2}$ vertical column density - ' + str(current_volcano.time_begin) + 'T' + str(current_volcano.hour) + ':' + str(current_volcano.minute) + ':' + str(current_volcano.second) + 'Z - ' + str(current_volcano.name) + ' - ' + str(method)
        textstr = "Total filled $SO_{2}$ mass of center volcano (" + current_volcano.name + "): " + "{:.2f}".format(mass_filled.sum() / 1000000) + " tons"
        plt.gcf().text(0.13, 0.06, textstr, fontsize=14)
        plt.title(title)

    if PRINT_MODE != 1:
        cba = plt.colorbar(cs)
        cba.set_label('$SO_{2}$ mole/$m^{2}$')

    if sem_inst == 0:
        cbb = plt.colorbar(cmask_filled, ticks=[0, 1])
        cbb.set_label('Filled mask')
    else:
        cbb = plt.colorbar(cmask_filled, ticks=FixedLocator(colorbar_labels))
        cbb.set_label('Associated volcano ID-s')

    cbc = plt.colorbar(cmask, ticks=[0, 1])
    cbc.set_label('$SO_{2}$ detection mask')

    # plt.show()

    return fig


def plot_ground_truth(current_volcano):
    """Plot the ground truth for a product. Takes in a Volcano class, returns the plot figure"""

    # OPEN FILE
    print(current_volcano.filename)
    fh = netCDF4.Dataset(current_volcano.filename, 'r', format="NETCDF4")

    PRINT_MODE = 0  # set to 1 if plotting for document, set to 0 if viewed on PC
    zoom_span = 6

    # PLOT
    # Set plot DPI high for Spyder plot
    if PRINT_MODE == 1:
        plt.rcParams['figure.dpi'] = 500
        fig = plt.figure(figsize=(4, 2))
    else:
        plt.rcParams['figure.dpi'] = 200
        fig = plt.figure(figsize=(16, 9))

    # GET DATA
    # - lat/lon
    lons = fh.groups['PRODUCT'].variables['longitude'][:][0, :, :]
    lats = fh.groups['PRODUCT'].variables['latitude'][:][0, :, :]
    # - other
    acqdate = fh.time_reference
    viewing_azimuth_angle = (
    fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['GEOLOCATIONS'].variables['viewing_azimuth_angle'][0, :, :])
    surface_altitude = (
    fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA'].variables['surface_altitude'][0, :, :])
    cloud_fraction = (
    fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA'].variables['cloud_fraction_crb'][0, :, :])
    cloud_height = (
    fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['INPUT_DATA'].variables['cloud_height_crb'][0, :, :])

    # PLOT
    m = Basemap(llcrnrlon=current_volcano.longitude - zoom_span / 2, llcrnrlat=current_volcano.latitude - zoom_span / 2,
                urcrnrlon=current_volcano.longitude + zoom_span / 2, urcrnrlat=current_volcano.latitude + zoom_span / 2,
                resolution='l', projection='cyl', lat_ts=40, lat_0=current_volcano.latitude,
                lon_0=current_volcano.longitude)

    gt_filename = ASSETS_DIRECTORY + "/ground_truth/" + str(os.path.basename(current_volcano.filename)) + ".csv"
    a = np.zeros(shape=lons.shape)
    if os.path.exists(gt_filename):
        a = np.genfromtxt(gt_filename, delimiter=',')

    u, ind = np.unique(a, return_inverse=True)
    b = ind.reshape((a.shape))

    labels = []
    print(np.unique(a))
    for element in np.unique(a):
        if not np.isnan(element):
            current_id = int(element)
        else:
            current_id = -1
        # print(current_id)
        if current_id == -1:
            current_label = str(current_id) + " - Unidentified"
            # print(current_label)
        elif current_id == 0:
            current_label = str(current_id) + " - No SO2"
            # print(current_label)
        else:
            current_label = str(current_id) + " - " + str(
                volcanolist_df.loc[volcanolist_df['id'] == current_id]['fullname'].values[0])
            # print(current_label)

        labels.append(current_label)

    # - SO2 detection mask
    m.pcolormesh(lons, lats, b, latlon=True, cmap=plt.cm.get_cmap('viridis', len(labels)), alpha=.5, shading='auto')

    for index, row in volcanolist_df.iterrows():
        plt.text(row['lon'], row['lat'], str(row['fullname']) + " (" + str(row['id']) + ")", clip_on=True, fontsize=6)
        plt.plot(row['lon'], row['lat'], '^b')

    m.drawparallels(np.arange(-85., 85., 2.), labels=[1, 0, 0, 0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 2.), labels=[0, 0, 0, 1], fontsize=10)
    try:
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
    except:
        print("No basemap lines drawn")

    cb = plt.colorbar(ticks=np.arange(len(u)))
    cb.ax.set_yticklabels(labels)
    if PRINT_MODE == 0:
        plt.title(current_volcano.time_begin + " - " + gt_filename)

    # plt.show()

    return fig


def plot_dbscan_clusters(current_volcano):
    """Plot the result of DBSCAN cluster generation. Takes in a Volcano class, returns the plot figure"""

    print(current_volcano.filename)
    # OPEN FILE
    fh = netCDF4.Dataset(current_volcano.filename, 'r', format="NETCDF4")

    PRINT_MODE = 0 # set to 1 if plotting for document, set to 0 if viewed on PC
    zoom_span = 6

    # PLOT
    # Set plot DPI high for Spyder plot
    if PRINT_MODE == 1:
        plt.rcParams['figure.dpi'] = 500
        fig = plt.figure(figsize=(4, 2))
    else:
        plt.rcParams['figure.dpi'] = 200
        fig = plt.figure(figsize=(16, 9))

    # GET DATA
    # - lat/lon
    lons = fh.groups['PRODUCT'].variables['longitude'][:][0, :, :]
    lats = fh.groups['PRODUCT'].variables['latitude'][:][0, :, :]

    # - SO2 maps
    SO2_PBL = fh.groups['PRODUCT'].variables['sulfurdioxide_total_vertical_column'][0, :, :]

    # - SO2 detection flag: 'sulfurdioxide_detection_flag'
    flag = fh.groups['PRODUCT'].groups['SUPPORT_DATA'].groups['DETAILED_RESULTS'].variables['sulfurdioxide_detection_flag'][0, :, :]
    mask = np.where(flag > 0, 1, 0)  # >> assign 1 for flag values above 0

    mask_filled = np.zeros(mask.shape)
    if mask.sum() >= 4:
        mask_filled = segmentation_algorithms.get_dbscan_label_map(SO2_PBL, mask)

    m = Basemap(llcrnrlon=current_volcano.longitude - zoom_span / 2,
                llcrnrlat=current_volcano.latitude - zoom_span / 2,
                urcrnrlon=current_volcano.longitude + zoom_span / 2,
                urcrnrlat=current_volcano.latitude + zoom_span / 2, resolution='l', projection='cyl',
                lat_ts=40, lat_0=current_volcano.latitude, lon_0=current_volcano.longitude)

    # - SO2 map
    vmin, vmax = 0, 0.01
    vcd2plot, vcd2plot_str = SO2_PBL, 'SO2_PBL'
    cs = m.pcolormesh(lons, lats, vcd2plot, latlon=True, vmin=vmin, vmax=vmax, cmap='CMRmap_r', alpha=.2, shading='auto')

    # Crop the clusters that are not on the plot, to eliminate non-present colorbar labels
    for iterator_row in range(len(mask_filled)):
        for iterator_column in range(len(mask_filled[0])):
            if lons[iterator_row, iterator_column] > current_volcano.longitude + zoom_span / 2 or lons[
                iterator_row, iterator_column] < current_volcano.longitude - zoom_span / 2:
                mask_filled[iterator_row, iterator_column] = 0
            if lats[iterator_row, iterator_column] > current_volcano.latitude + zoom_span / 2 or lats[
                iterator_row, iterator_column] < current_volcano.latitude - zoom_span / 2:
                mask_filled[iterator_row, iterator_column] = 0

    # Count the remaining unique elements after cropping the image
    mask_filled = mask_filled.astype(int)
    number_of_remained_unique_elements_after_crop = len(np.unique(mask_filled))
    mask_filled = mask_filled.astype(np.float32)

    colorbar_labels = []

    for element in np.unique(mask_filled[~np.isnan(mask_filled)]):
        current_id = int(element)
        if current_id != 0:
            current_label = int(current_id)
            colorbar_labels.append(current_label)
    if len(colorbar_labels) == 0:
        colorbar_labels.append(0)
        mycmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["w", "w"])
        cmask_filled = m.pcolormesh(lons, lats, mask_filled, latlon=True, cmap=mycmap, alpha=.8, shading='auto')
    else:
        mask_filled[mask_filled == 0] = np.nan  # to make the zeros transparent in pcolormesh
        mycmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["C0", "yellow", "C1", "C3", "C2"], number_of_remained_unique_elements_after_crop - 1)
        cmask_filled = m.pcolormesh(lons, lats, mask_filled, latlon=True, cmap=mycmap, alpha=.8, shading='auto')

    # - SO2 detection mask
    cmask = m.pcolormesh(lons, lats, mask.astype(np.float32), latlon=True, vmin=0, vmax=1, cmap=plt.cm.get_cmap('Greys', 2), alpha=.1, shading='auto')

    cbb = plt.colorbar(cmask_filled, ticks=FixedLocator(colorbar_labels))
    cbb.set_label('Cluster labels')

    for other_volcanoindex, other_row in volcanolist_df.iterrows():
        plt.text(other_row['lon'], other_row['lat'], str(other_row['fullname']) + " (" + str(other_row['id']) + ")", clip_on=True, fontsize=6)
        plt.plot(other_row['lon'], other_row['lat'], '^b')

    m.drawparallels(np.arange(-85., 85., 2.), labels=[1, 0, 0, 0], fontsize=10)
    m.drawmeridians(np.arange(-180., 181., 2.), labels=[0, 0, 0, 1], fontsize=10)
    try:
        m.drawcoastlines()
        m.drawstates()
        m.drawcountries()
    except:
        print("No basemap lines drawn")

    return fig
