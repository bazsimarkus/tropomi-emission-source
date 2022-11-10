# -----------------------------------------------------------------------------------------------
# SET THIS TO 1 IF AN EXTERNAL HARD DRIVE CALLED "Main Volume" IS MOUNTED
# SET THIS TO 0. IF THE ASSETS ARE PLACED IN THE REPOSITORY ROOT IN A FOLDER NAMED "assets"

EXTERNAL_DRIVE = 0

# -----------------------------------------------------------------------------------------------

# DON'T EDIT AFTER THIS

molar_mass_SO2 = 64.0638
pixel_area_m2 = 19250000

import getpass
import os

if EXTERNAL_DRIVE == 1:
    username = str(getpass.getuser())
    ASSETS_DIRECTORY = "/media/" + username + "/Main Volume"
else:
    dirname = os.path.dirname(__file__)
    ASSETS_DIRECTORY = os.path.join(dirname, 'assets')
    print(ASSETS_DIRECTORY)