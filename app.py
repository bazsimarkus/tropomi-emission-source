# -----------------------------------------------------------
# app.py
# Contains the API to run and test the algorithms. Offers an alternative and intuitive running method to main.py
# -----------------------------------------------------------

import io
import warnings

import cv2
import numpy as np
from fastapi import FastAPI, Query
from starlette.responses import StreamingResponse

import plotting_functions
import segmentation_algorithms
from Volcano import *
from calculate_metrics import get_metrics_response

warnings.filterwarnings("ignore")

app = FastAPI(title='TROPOMI SO2 API', description='Automatic Retrieval of Volcanic SO2 Emission Source from TROPOMI Products')

volcanolist_df = pd.read_csv('volcano_coordinates.csv')
volcano_name_list = sorted(volcanolist_df['name'].tolist(), key=str.lower)
method_list = ["radius_100km", "floodfill", "floodfill_withwind", "dbscan", "dbscan_withwind", "reverse_trajectory_dbscan", "multi_dbscan", "multi_reverse_trajectory_dbscan"]


@app.get("/return_so2_mass_of_volcano", tags=["Return and plot SO2 mass of a volcano"], description="Return the associated SO2 mass of one volcano - JSON response", responses={"200":{"description":"Successful Response","content":{"application/json":{"example":{"data":[{"name":"sangay","fullname":"Sangay","datetime":"2021-08-15T18:41:35Z","latitude":-2.005,"longitude":-78.341,"elevation":5286,"method":"floodfill","id":55,"filename":"S5P_NRTI_L2__SO2____20210815T184135_20210815T184635_19896_02_020201_20210815T193716.nc","so2_mass_value":537.8432463777182,"so2_mass_unit":"tons","ground_truth_available":1,"intersection_over_union":0.8111782477341389,"oversegmentation":0,"undersegmentation":0.18882175226586106,"area_fit_index":-0.23277467411545624,"root_mean_square":0.1335171414627167,"true_positive_pixels":537,"true_negative_pixels":212,"false_positive_pixels":0,"false_negative_pixels":125,"precision":1,"recall":0.8111782477341389,"f1_score":0.8957464553794828,"accuracy":0.8569794050343249,"specificity":1,"ground_truth_pixels_for_current_volcano":662,"so2_pixels_on_image":874}]}}}}})
async def return_so2_mass_of_volcano(volcano_name: str = Query(enum=volcano_name_list, description="The lowercase name of the volcano, that is present in the volcano_coordinates.csv list. Example: **sangay**"), volcano_date: str = Query(regex="[0-9]{4}-[0-9]{2}-[0-9]{2}", description="The date of the satellite image, that we want to run the segmentation on. Example: **2021-08-15**"), method: str = Query(enum=method_list, description="The method that we want to use for the segmentation. Example: **dbscan**")):
    """Return the associated SO2 mass and the resulting metrics for one volcano - JSON response"""

    current_volcano = Volcano(volcano_name, volcano_date)

    function_to_use = getattr(segmentation_algorithms, 'run_' + str(method))
    result = function_to_use(current_volcano)

    response = get_metrics_response(*result)

    return response


@app.get("/plot_so2_mass_of_volcano", tags=["Return and plot SO2 mass of a volcano"], description="Return the plot of the associated SO2 plumes for one volcano as a .PNG image", response_class=StreamingResponse)
async def plot_so2_mass_of_volcano(volcano_name: str = Query(enum=volcano_name_list, description="The lowercase name of the volcano, that is present in the volcano_coordinates.csv list. Example: **sangay**"), volcano_date: str = Query(regex="[0-9]{4}-[0-9]{2}-[0-9]{2}", description="The date of the satellite image, that we want to run the segmentation on. Example: **2021-08-15**"), method: str = Query(enum=method_list, description="The method that we want to use for the segmentation. Example: **dbscan**")):
    """Return the plot of the associated SO2 mass for one volcano - Image response"""

    current_volcano = Volcano(volcano_name, volcano_date)

    function_to_use = getattr(segmentation_algorithms, 'run_' + str(method))
    result = function_to_use(current_volcano)

    fig = plotting_functions.plot_segmentation_result(*result)
    fig.canvas.draw()
    # convert canvas to image
    img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    # img is rgb, convert to opencv's default bgr
    img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
    res, im_png = cv2.imencode(".png", img)
    return StreamingResponse(io.BytesIO(im_png.tobytes()), media_type="image/png")


@app.get("/plot_ground_truth", tags=["Utilities"], description="This endpoint plots the ground truth for the image. If the ground truth doesn't exist, it assumes that there is an empty ground truth.")
async def plot_ground_truth(volcano_name: str = Query(enum=volcano_name_list, description="The lowercase name of the volcano, that is present in the volcano_coordinates.csv list. Example: **sangay**"), volcano_date: str = Query(regex="[0-9]{4}-[0-9]{2}-[0-9]{2}", description="The date of the satellite image, that we want to run the segmentation on. Example: **2021-08-15**")):
    """Return the plot of the ground truth for one volcano - Image response"""

    current_volcano = Volcano(volcano_name, volcano_date)

    function_to_use = getattr(plotting_functions, 'plot_ground_truth')
    result = function_to_use(current_volcano)

    fig = result
    fig.canvas.draw()
    # convert canvas to image
    img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    # img is rgb, convert to opencv's default bgr
    img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
    res, im_png = cv2.imencode(".png", img)
    return StreamingResponse(io.BytesIO(im_png.tobytes()), media_type="image/png")


@app.get("/plot_dbscan_clusters", tags=["Utilities"], description="This endpoint plots the generated DBSCAN clusters with their labels for the image.")
async def plot_dbscan_clusters(volcano_name: str = Query(enum=volcano_name_list, description="The lowercase name of the volcano, that is present in the volcano_coordinates.csv list. Example: **sangay**"), volcano_date: str = Query(regex="[0-9]{4}-[0-9]{2}-[0-9]{2}", description="The date of the satellite image, that we want to run the segmentation on. Example: **2021-08-15**")):
    """Return the plot of generated DBSCAN clusters for one volcano - Image response"""

    current_volcano = Volcano(volcano_name, volcano_date)

    function_to_use = getattr(plotting_functions, 'plot_dbscan_clusters')
    result = function_to_use(current_volcano)

    fig = result
    fig.canvas.draw()
    # convert canvas to image
    img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    # img is rgb, convert to opencv's default bgr
    img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
    res, im_png = cv2.imencode(".png", img)
    return StreamingResponse(io.BytesIO(im_png.tobytes()), media_type="image/png")
