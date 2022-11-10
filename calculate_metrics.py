# -----------------------------------------------------------
# calculate_metrics.py
# Contains the function for relevant metric calculation to evaluate the algorithms
# -----------------------------------------------------------

import math

import numpy as np

from CONFIG_PARAMETERS import molar_mass_SO2, pixel_area_m2
from Volcano import *


def get_metrics_response(sem_inst, mask_filled, gt_data, current_volcano, lons, lats, mask, SO2_PBL, mass_filled, method, trajectories=None):
    """Calculate relevant metrics to evaluate the algorithms. Takes in the result of an algorithm function in segmentation_algorithms.py, returns an object containing the results"""

    data = {"data": []}

    if sem_inst == 0:
        volcanoes_to_check = [current_volcano.id]
    else:
        volcanoes_to_check = [current_volcano.id]
        # volcanoes_to_check = np.unique(mask_filled)

    for current_volcano_id in volcanoes_to_check:
        if current_volcano_id != 0:
            current_volcano_data = volcanolist_df.loc[volcanolist_df['id'] == int(current_volcano_id)]

            current_mask_filled = np.where(mask_filled == current_volcano_id, current_volcano_id, 0)
            current_mask_filled_binary = np.where(mask_filled == current_volcano_id, 1, 0)
            center_volcano_mask = np.where(mask_filled == current_volcano_id, 1, 0)
            mass_filled = np.multiply(SO2_PBL, center_volcano_mask)
            mass_filled = mass_filled * molar_mass_SO2 * pixel_area_m2
            so2_mass = mass_filled.sum() / 1000000

            current_gt_data = gt_data
            so2_pixels_on_image = int(np.count_nonzero(mask == 1))

            if current_gt_data is not None:
                # Mask GT for current volcano
                current_gt_data = np.where(current_gt_data == current_volcano_id, 1, 0)  # >> assign 1 for flag values above 0
                ground_truth_available = 1

                true_positive_pixels = int(np.sum(np.logical_and(current_mask_filled_binary == 1, current_gt_data == 1)))
                true_negative_pixels = int(np.sum(np.logical_and(current_mask_filled_binary == 0, current_gt_data == 0))) - int(np.count_nonzero(mask == 0))
                false_positive_pixels = int(np.sum(np.logical_and(current_mask_filled_binary == 1, current_gt_data == 0)))
                false_negative_pixels = int(np.sum(np.logical_and(current_mask_filled_binary == 0, current_gt_data == 1)))

                ground_truth_pixels_for_current_volcano = int(current_gt_data.sum())

                # CALCULATION OF SEGMENTATION METRICS
                # Intersection over Union

                if (true_positive_pixels + false_positive_pixels + false_negative_pixels) != 0:
                    intersection_over_union = true_positive_pixels / (true_positive_pixels + false_positive_pixels + false_negative_pixels)
                else:
                    intersection_over_union = None

                # Area Fit Index

                if (true_positive_pixels + false_positive_pixels) != 0:
                    area_fit_index = (false_positive_pixels - false_negative_pixels) / (true_positive_pixels + false_positive_pixels)
                else:
                    area_fit_index = None

                # Oversegmentation, undersegmentation

                if (true_positive_pixels + false_positive_pixels) != 0:
                    oversegmentation = 1 - true_positive_pixels / (true_positive_pixels + false_positive_pixels)
                else:
                    oversegmentation = None

                if (true_positive_pixels + false_negative_pixels) != 0:
                    undersegmentation = 1 - true_positive_pixels / (true_positive_pixels + false_negative_pixels)
                else:
                    undersegmentation = None

                # Root Mean Square

                if oversegmentation is not None and undersegmentation is not None:
                    root_mean_square = math.sqrt((oversegmentation * oversegmentation + undersegmentation * undersegmentation) / 2)
                else:
                    root_mean_square = None

                # Precision, recall, F1, accuracy, specificity

                if (true_positive_pixels + false_positive_pixels) != 0:
                    precision = true_positive_pixels / (true_positive_pixels + false_positive_pixels)
                else:
                    precision = None

                if (true_positive_pixels + false_negative_pixels) != 0:
                    recall = true_positive_pixels / (true_positive_pixels + false_negative_pixels)
                else:
                    recall = None

                if precision is not None and recall is not None:
                    f1_score = (2 * precision * recall) / (precision + recall)
                else:
                    f1_score = None

                if (true_positive_pixels + true_negative_pixels + false_positive_pixels + false_negative_pixels) != 0:
                    accuracy = (true_positive_pixels + true_negative_pixels) / (true_positive_pixels + true_negative_pixels + false_positive_pixels + false_negative_pixels)
                else:
                    accuracy = None

                if (true_negative_pixels + false_positive_pixels) != 0:
                    specificity = true_negative_pixels / (true_negative_pixels + false_positive_pixels)
                else:
                    specificity = None

            else:
                ground_truth_available = 0
                intersection_over_union = None
                oversegmentation = None
                undersegmentation = None
                area_fit_index = None
                root_mean_square = None

                true_positive_pixels = None
                true_negative_pixels = None
                false_positive_pixels = None
                false_negative_pixels = None

                ground_truth_pixels_for_current_volcano = None

                precision = None
                recall = None
                f1_score = None
                accuracy = None
                specificity = None

            response = {
                    "name": str(current_volcano_data['name'].to_string(index=False)),
                    "fullname": str(current_volcano_data['fullname'].to_string(index=False)),
                    "datetime": str(current_volcano.time_begin) + "T" + str(current_volcano.hour) + ":" + str(current_volcano.minute) + ":" + str(current_volcano.second) + "Z",
                    "latitude": float(current_volcano_data['lat']),
                    "longitude": float(current_volcano_data['lon']),
                    "elevation": int(current_volcano_data['elevation']),
                    "method": method,
                    "id": int(current_volcano_id),
                    "filename": os.path.basename(current_volcano.filename),
                    "so2_mass_value": so2_mass,
                    "so2_mass_unit": "tons",
                    "ground_truth_available": ground_truth_available,
                    "intersection_over_union": intersection_over_union,
                    "oversegmentation": oversegmentation,
                    "undersegmentation": undersegmentation,
                    "area_fit_index": area_fit_index,
                    "root_mean_square": root_mean_square,
                    "true_positive_pixels": true_positive_pixels,
                    "true_negative_pixels": true_negative_pixels,
                    "false_positive_pixels": false_positive_pixels,
                    "false_negative_pixels": false_negative_pixels,
                    "precision": precision,
                    "recall": recall,
                    "f1_score": f1_score,
                    "accuracy": accuracy,
                    "specificity": specificity,
                    "ground_truth_pixels_for_current_volcano": ground_truth_pixels_for_current_volcano,
                    "so2_pixels_on_image": so2_pixels_on_image
                }

            data['data'].append(response)

    return data
