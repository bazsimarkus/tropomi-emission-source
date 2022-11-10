# -----------------------------------------------------------
# main.py
# Runs the result generation using the predefined parameters in the CONFIG files
# -----------------------------------------------------------

import csv
import datetime
import gc
import json
import os
import warnings

import segmentation_algorithms
from Volcano import *
from calculate_metrics import get_metrics_response
from plotting_functions import plot_segmentation_result

warnings.filterwarnings("ignore")

# Datetime object for the distinct log file name
now = datetime.datetime.now()
dt_string = now.strftime("%Y%m%dT%H%M%S")

volcano_list_df = pd.read_csv('CONFIG_INPUT_LIST.csv')

OUTPUT_CSV_OR_JSON = 1
fnames = ['name','fullname','id','datetime','latitude','longitude','elevation','method','filename','so2_mass_value','so2_mass_unit','so2_pixels_on_image','ground_truth_available','ground_truth_pixels_for_current_volcano','true_positive_pixels','false_positive_pixels','true_negative_pixels','false_negative_pixels','intersection_over_union','oversegmentation','undersegmentation','area_fit_index','root_mean_square','accuracy','precision','recall','f1_score','specificity']

if OUTPUT_CSV_OR_JSON == 1:
    output_filename = 'output/result_' + dt_string + '.csv'
    with open(output_filename, 'a') as filedata:
        writer = csv.writer(filedata, delimiter=',')
        writer.writerow(fnames)


for index, row in volcano_list_df.iterrows():
    try:
        current_volcano = Volcano(row['name'], row['date'])

        function_to_use = getattr(segmentation_algorithms, 'run_' + str(row['method']))

        result = function_to_use(current_volcano)

        response = get_metrics_response(*result)

        if OUTPUT_CSV_OR_JSON == 0:
            output_filename = 'output/result_' + dt_string + '.json'
            if not os.path.exists(output_filename):
                data = {"data": []}
                with open(output_filename, 'w') as outputfile:
                    json.dump(data, outputfile)

            with open(output_filename, 'r+') as outputfile:
                # First we load existing data into a dict.
                file_data = json.load(outputfile)
                # Join new_data with file_data inside emp_details
                file_data["data"].append(response['data'])
                # Sets file's current position at offset.
                outputfile.seek(0)
                # convert back to json.
                json.dump(file_data, outputfile, indent=4)
        else:
            output_filename = 'output/result_' + dt_string + '.csv'
            with open(output_filename, 'a') as filedata:
                writer = csv.DictWriter(filedata, delimiter=',', fieldnames=fnames)
                writer.writerow(response['data'][0])

        fig = plot_segmentation_result(*result)
        fig.savefig('plots/' + str(row['name']) + '_' + str(row['date']) + '_' + str(row['method']) + '_' + dt_string + '.png', bbox_inches='tight')

        gc.collect()
    except Exception as e:
        print("Exception: " + str(e))
