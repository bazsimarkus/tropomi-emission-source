# -----------------------------------------------------------
# trajectory_algorithms.py
# Contains all functions to generate forward and reverse trajectories for the wind-based algorithms in segmentation_algorithms.py
# -----------------------------------------------------------

import pysplit
import os
from CONFIG_PARAMETERS import ASSETS_DIRECTORY


def get_forward_trajectory(current_volcano):
    """Generate a forward trajectory. Takes in a Volcano class, returns the lat/lon points of the trajectory"""

    basename = 'volcano'

    working_dir = ASSETS_DIRECTORY + '/hysplit/working'
    storage_dir = ASSETS_DIRECTORY + '/trajectories/' + basename + "/"
    meteo_dir = ASSETS_DIRECTORY + '/gdas'

    years = [current_volcano.year]
    months = [current_volcano.month]
    day = current_volcano.day
    hours = [current_volcano.hour - 1]
    altitudes = [current_volcano.elevation]
    location = (current_volcano.latitude, current_volcano.longitude)

    runtime = 12

    pysplit.generate_bulktraj(basename, working_dir, storage_dir, meteo_dir,
                              years, months, hours, altitudes, location, runtime,
                              monthslice=slice(day - 1, day, 1), get_reverse=False,
                              get_clipped=True)
    list_of_files = [storage_dir + name for name in os.listdir(storage_dir) if
                     os.path.isfile(os.path.join(storage_dir, name))]

    latest_file = max(list_of_files, key=os.path.getmtime)
    hysplit_filename = str(os.path.basename(latest_file))

    # print(hysplit_filename)

    test_data = pysplit.load_hysplitfile(storage_dir + hysplit_filename)
    xtraj = []
    ytraj = []

    for current_point in test_data[1]:
        lat_coord = current_point[0]
        lon_coord = current_point[1]
        xtraj.append(lat_coord)
        ytraj.append(lon_coord)

    return xtraj, ytraj


def get_reverse_trajectory(latitude, longitude, year, month, day, hour, elevation):
    """Generate a reverse trajectory. Takes in latitude, longitude, year, month, day, hour, and elevation to start from, returns the lat/lon points of the trajectory"""

    basename = 'volcano'

    working_dir = ASSETS_DIRECTORY + '/hysplit/working'
    storage_dir = ASSETS_DIRECTORY + '/trajectories/' + basename + "/"
    meteo_dir = ASSETS_DIRECTORY + '/gdas'

    years = [year]
    months = [month]
    day = day
    hours = [hour]
    altitudes = [elevation]
    location = (latitude, longitude)

    runtime = -12

    pysplit.generate_bulktraj(basename, working_dir, storage_dir, meteo_dir,
                              years, months, hours, altitudes, location, runtime,
                              monthslice=slice(day - 1, day, 1), get_reverse=False,
                              get_clipped=True)
    list_of_files = [storage_dir + name for name in os.listdir(storage_dir) if
                     os.path.isfile(os.path.join(storage_dir, name))]

    latest_file = max(list_of_files, key=os.path.getmtime)
    hysplit_filename = str(os.path.basename(latest_file))

    # print(hysplit_filename)

    test_data = pysplit.load_hysplitfile(storage_dir + hysplit_filename)
    xtraj = []
    ytraj = []

    for current_point in test_data[1]:
        lat_coord = current_point[0]
        lon_coord = current_point[1]
        xtraj.append(lat_coord)
        ytraj.append(lon_coord)

    return xtraj, ytraj
