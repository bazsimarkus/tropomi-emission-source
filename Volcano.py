# -----------------------------------------------------------
# Volcano.py
# Contains the Volcano class, which contains all relevant information to a volcanic scenario (name, date, location, filename)
# -----------------------------------------------------------

import os

import pandas as pd

from CONFIG_PARAMETERS import ASSETS_DIRECTORY

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

volcanolist_df = pd.read_csv('volcano_coordinates.csv')


class Volcano:
    """Volcano class implementation"""
    def __init__(self, name, time_begin):
        self.name = name
        self.time_begin = time_begin

        self.year = int(time_begin.split('-')[0])
        self.month = int(time_begin.split('-')[1])
        self.day = int(time_begin.split('-')[2])

        self.fullname = str(volcanolist_df.loc[volcanolist_df['name'] == name]['fullname'].to_string(index=False))
        self.latitude = float(volcanolist_df.loc[volcanolist_df['name'] == name]['lat'])
        self.longitude = float(volcanolist_df.loc[volcanolist_df['name'] == name]['lon'])
        self.elevation = int(volcanolist_df.loc[volcanolist_df['name'] == name]['elevation'])
        self.id = float(volcanolist_df.loc[volcanolist_df['name'] == name]['id'])

        # Connect to the db file
        tropomi_files = pd.read_csv('CONFIG_TROPOMI_file_archive.csv')

        # Get path to file in local database
        dataset_root = ASSETS_DIRECTORY
        print(dataset_root)
        self.filename = dataset_root + tropomi_files.loc[(tropomi_files['name'] == name) & (tropomi_files['date'] == time_begin)]['filepath'].to_string(index=False)
        self.hour = int(str(os.path.basename(self.filename)).replace("S5P_NRTI_L2__SO2____","")[9:11])
        self.minute = int(str(os.path.basename(self.filename)).replace("S5P_NRTI_L2__SO2____","")[11:13])
        self.second = int(str(os.path.basename(self.filename)).replace("S5P_NRTI_L2__SO2____","")[13:15])
