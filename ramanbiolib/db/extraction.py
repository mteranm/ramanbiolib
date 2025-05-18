"""
Functions to extract spectra traces from articles figures plots.
"""

import math
import pytesseract
import cv2 as cv
import numpy as np
import pandas as pd
import plotly.express as px

from scipy import interpolate
from scipy.signal import find_peaks, savgol_filter

def __get_y_axis_ticks(im, axis_contour):
    axis_min_y_value = axis_contour[:,:,1].max()
    ticks_str = pytesseract.image_to_string(im[axis_min_y_value:, :, :])
    print(f"""Axis extracted labels (OCR):
{ticks_str}
    """)
    lines = ticks_str.split('\n')
    lines = list(filter(lambda x: len(x.strip()) > 0, lines))
    return [ int(value) for value in lines[0].split(' ')]

def __create_wavenumbers_vector_simple(raman_start, step, vector_len):
    raman_shift = []
    for i in range(vector_len):
        raman_shift.append(raman_start + i*step)
    return raman_shift

def __create_wavenumbers_vector(raman_pixels, raman_ticks, trace_min, trace_max):
    wavenumbers = []
    pixels = []

    raman_pixels = list(filter(lambda x: x < trace_max, raman_pixels))

    for i in range(len(raman_pixels)):
        if i+1 < len(raman_pixels):
            step = (raman_ticks[i+1]-raman_ticks[i])/(raman_pixels[i+1]-raman_pixels[i])
            if i == 0:
                # First trace part before first tick
                vector_len =  raman_pixels[0] - trace_min
                raman_ticks[i]-raman_ticks[i+1]
                raman_start = raman_ticks[0] - vector_len*step
                wavenumbers.extend(
                    __create_wavenumbers_vector_simple(raman_start, step, vector_len)
                )
            # Between ticks
            wavenumbers.extend(
                __create_wavenumbers_vector_simple(
                    raman_ticks[i], 
                    step, 
                    raman_pixels[i+1]-raman_pixels[i]
                )
            )
        else:
            # Last trace part after last tick
            vector_len =   trace_max - raman_pixels[i] + 1
            step = (raman_ticks[i]-raman_ticks[i-1])/(raman_pixels[i]-raman_pixels[i-1])
            raman_start = raman_ticks[i]
            wavenumbers.extend(
                __create_wavenumbers_vector_simple(raman_start, step, vector_len)
            )
    return wavenumbers

def __extract_wavenumbers(trace, axis_contour, raman_ticks, ticks_down, prominence_ticks, 
    ignore_first_tick, im):
    # Group and reduce axis points
    points = { v:[] for v in np.unique(axis_contour.reshape((axis_contour.shape[0], 2))[:,0])}
    for point in axis_contour.reshape((axis_contour.shape[0], 2)):
        points[point[0]].append(point[1])
    min_axes_values = np.array([ 
        [x, np.max(y) if ticks_down else np.min(y)] 
        for x, y in points.items()
    ])
    # Find peaks to find aixs ticks positions
    aixs_line = min_axes_values[:,1] if ticks_down else  min_axes_values[:,1]*(-1)

    px.line(aixs_line).write_image("logs/axis-peaks-line.png")

    if ignore_first_tick:
        # Somtimes first tick is the y axis so we need to remove it
        peaks, _ = find_peaks(aixs_line)
        if peaks[0] < 100: # Sometimes the first peak can be detected or not depending on the run (starting from 0 or from the peak)
            ignoring_shift = peaks[0]+10
            aixs_line = aixs_line[ignoring_shift:]
            px.line(aixs_line).write_image("logs/axis-peaks-line2.png")
        else:
            ignoring_shift = peaks[0]-10
            aixs_line = aixs_line[ignoring_shift:]
            px.line(aixs_line).write_image("logs/axis-peaks-line2.png")

    aixs_norm =np.pad(
        (aixs_line - np.min(aixs_line))/
        (np.max(aixs_line) - np.min(aixs_line)),
         (10, 10)
    )
    px.line(aixs_norm).write_image("logs/axis-peaks.png")

    peaks, _ = find_peaks(aixs_norm, 
        distance=20, prominence=prominence_ticks
    )

    if ignore_first_tick:
        # we need to correct for first tick ignore when finidng peaks
        print(peaks)
        peaks = [ p + ignoring_shift for p in peaks]
        print(peaks)

    raman_pixels = np.pad(min_axes_values[:,0], (10, 10))[peaks]

    for raman_pixel in raman_pixels:
        cv.line(im, (raman_pixel, 0), (raman_pixel, im.shape[0]), (255, 0, 0) , 2) 
    px.imshow(im).write_image("logs/axis-ticks.png")
    
    print(f"Wavenumbers\n ticks: {raman_ticks}\n pixels: {raman_pixels}")

    if len(raman_pixels) != len(raman_ticks):
        raise Exception("Wrong ticks pixels extraction")

    # Create wavenumber vector from ticks positions
    return __create_wavenumbers_vector(raman_pixels, raman_ticks, trace.min(), trace.max())
 

def __extract_trace(
    img_shape, contour, subfigure, component_name, axis_contour, axis_y_ticks,
    ticks_down, prominence_ticks, ignore_first_tick, im
):
    # Group and reduce trace points
    points = {v: [] for v in np.unique(contour.reshape((contour.shape[0], 2))[:, 0])}
    for point in contour.reshape((contour.shape[0], 2)):
        points[point[0]].append(point[1])

    # Create dataframe row
    trace = pd.DataFrame(
        [
            {
                "x": x,
                "intensity": (np.min(y) - img_shape) * -1,
                "subfigure": subfigure,
                "component_name": component_name,
            }
            for x, y in points.items()
        ]
    ).sort_values(by="x")

    # Extract the wavenumbers for each intensity point
    x = np.array(__extract_wavenumbers(trace["x"], axis_contour, axis_y_ticks, 
        ticks_down, prominence_ticks, ignore_first_tick, im)
    )
    y = trace["intensity"].to_numpy()
    print(f"x len: {len(x)}, y len: {len(y)}")
    interp_trace = pd.DataFrame(
        {
            "wavenumbers": x,
            "intensity": y,
            "component_name": [component_name for v in range(len(x))],
        }
    )
    return interp_trace


def __get_contours_info(imgray, threshold):
    ret, thresh = cv.threshold(imgray, threshold, 255, cv.THRESH_BINARY_INV)
    contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)

    contours_info = []
    for i, c in enumerate(contours):
        max_point = c.reshape(c.shape[0], 2)[:, 1].max()
        contours_info.append({"size": c.shape[0], "max_point": max_point, "contour": c})
    contours_info.sort(key=lambda x: x["size"], reverse=True)
    return contours_info

def __get_axis_max_point(imgray, img_copy):
    max_point = 0
    middle_point = int(imgray.shape[0]/2)
    lines = cv.HoughLinesP(cv.Canny(imgray[middle_point:,:], 80, 120),
        rho = 1,theta = 1*np.pi/4, 
        threshold = 100,
        minLineLength = 300,
        maxLineGap = 10
    )
    for i, line in enumerate(lines):
        line = line[0]
        pt1 = (line[0],line[1]+middle_point+max_point)
        pt2 = (line[2],line[3]+middle_point+max_point)
        img_copy = img_copy.copy()
        cv.line(img_copy, pt1, pt2, (0,0,255), 3)
        max_point = np.min([line[1], line[3]])
        print(line[2]-line[0])
        if line[2]-line[0] > 0:
            print(line)
            cv.imwrite("logs/axis_line.png", img_copy)
            break
    return middle_point+max_point

def extract_traces_and_axis(figure_info, threshold=127, axis_y_ticks=None,
        same_contour_axis_trace=False, ticks_down=True, prominence_ticks=0.5,
        ignore_first_tick=False, horizontal_crop=0):
    """Extract traces from a figure plot, using opencv contour detection after binary thresholding.

    Parameters
    ----------
    figure_info : dict
        The figure information dict in the format: { "image_path": "{path}", "subfigures": [ { "text": "{bottom-trace-name}"}, ...]}
    threshold : int
        The binary threshold value
    axis_y_ticks : list
        The list of the ticks wavenumbers values. If None OCR is tried for extract the values from the figure.
    same_contour_axis_trace : bool
        A flag to indicate if the bottom trace and axis are overlapping.
    ticks_down : bool
        A flag to indicate if the x-axis ticks are down the axis.
    prominence_ticks : bool
        The minimum prominence used in the x-axis ticks peaks detection. Useful to fine-tune the ticks detection.
    ignore_first_tick : bool
        A flag to indicate if the first tick detected needs to be avoided. Useful in cases where the ticks are over the x-axis and therefore the y-axis can be confused with a tick.        
    horizontal_crop : int
        The horizontal crop min value in pixels. Default 0, no crop. Useful when is necessary to do not consider some region of the plot image.
    """
    traces = []
    subfigures_num = len(figure_info["subfigures"])

    # Load image and find contours
    im = cv.imread(figure_info["image_path"])[:,horizontal_crop:,:]
    imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)

    # Get main contours
    contours_info = __get_contours_info(imgray, threshold)
    contours_info.sort(key=lambda x: x["size"], reverse=True)
    traces_axis_contour = contours_info[: subfigures_num + 1]
    traces_axis_contour.sort(key=lambda x: x["max_point"], reverse=True)

    print(f"Found {len(contours_info)} contours")

    # Get axis contour
    axis_contour = traces_axis_contour[0]
    print(f"Axis countour size: {axis_contour['size']}")
    px.imshow(cv.drawContours(im.copy(), [axis_contour['contour']], 0, (0,255,0), 3)).write_image("logs/axis.png")

    # Get axis ticks
    if axis_y_ticks is None:
        axis_y_ticks = __get_y_axis_ticks(im, axis_contour["contour"])

    if same_contour_axis_trace:
        axis_max_point = __get_axis_max_point(imgray, im.copy())
        print(axis_max_point)
        px.imshow(im[:axis_max_point,:,:]).write_image(f"logs/traces-noaxis.png")
        contours_info = __get_contours_info(imgray, 120)
        contours_info = __get_contours_info(imgray[:axis_max_point,:], threshold)
        contours_info.sort(key=lambda x: x["size"], reverse=True)
        traces_axis_contour = contours_info[: subfigures_num]
        traces_axis_contour.sort(key=lambda x: x["max_point"], reverse=True)

    # Get traces contour
    traces_contour = traces_axis_contour[:subfigures_num] if same_contour_axis_trace else traces_axis_contour[1:]
    for i, trace_c in enumerate(traces_contour):
        print(f"Trace {i} countour size: {trace_c['size']}")
        px.imshow(cv.drawContours(im.copy(), [trace_c['contour']], 0, (0,255,0), 3)).write_image(f"logs/trace{i}.png")

    # Create traces
    for i, trace_contour in enumerate(traces_contour):
        component_name = figure_info["subfigures"][i]["text"]
        print(f"Processing trace {component_name} ...")
        trace_df = __extract_trace(
            im.shape[0],
            trace_contour["contour"],
            'a',
            component_name,
            axis_contour["contour"],
            axis_y_ticks,
            ticks_down,
            prominence_ticks,
            ignore_first_tick,
            im.copy()
        )
        traces.append(trace_df)
    return traces

def extract_colored_traces_and_axis(figure_info, threshold=127, axis_y_ticks=None,
        same_contour_axis_trace=False, ticks_down=True, prominence_ticks=0.5,
        ignore_first_tick=False, horizontal_crop_min=0, horizontal_crop_max=-1):
    """Extract colored traces from a figure plot, using opencv contour detection after thresholding for color images.

    Parameters
    ----------
    figure_info : dict
        The figure information dict in the format: { "image_path": "{path}", "subfigures": [ { "text": "{bottom-trace-name}"}, ...]}
    threshold : int
        The binary threshold value
    axis_y_ticks : list
        The list of the ticks wavenumbers values. If None OCR is tried for extract the values from the figure.
    same_contour_axis_trace : bool
        A flag to indicate if the bottom trace and axis are overlapping.
    ticks_down : bool
        A flag to indicate if the x-axis ticks are down the axis.
    prominence_ticks : bool
        The minimum prominence used in the x-axis ticks peaks detection.
    ignore_first_tick : bool
        A flag to indicate if the first tick detected needs to be avoided. Useful in cases where the ticks are over the x-axis and therefore the y-axis can be confused with a tick.        
    horizontal_crop_min : int
        The horizontal crop min value in pixels. Default 0, no crop. Useful when is necessary to do not consider some region of the plot image.
    horizontal_crop_max : int
        The horizontal crop max value in pixels. Default -1, no crop. Useful when is necessary to do not consider some region of the plot image.
    """
    traces = []
    subfigures_num = len(figure_info["subfigures"])

    # Load image and find contours
    im = cv.imread(figure_info["image_path"])[:,horizontal_crop_min:horizontal_crop_max,:]
    imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)

    # Get axis contour
    contours_info = __get_contours_info(imgray, threshold)
    contours_info.sort(key=lambda x: x["size"], reverse=True)
    axis_contour = contours_info[0]
    print(f"Axis countour size: {axis_contour['size']}")
    px.imshow(cv.drawContours(im.copy(), [axis_contour['contour']], 0, (0,255,0), 3)).write_image("logs/axis.png")

    # Get axis ticks
    if axis_y_ticks is None:
        axis_y_ticks = __get_y_axis_ticks(im, axis_contour["contour"])

    # Get traces contour
    hsv = cv.cvtColor(im, cv.COLOR_RGB2HSV)
    lower_gray = np.array([0, 0, 0]) # Define range of gray color in HSV
    upper_gray = np.array([255, 10, 255])    
    imgray = cv.inRange(hsv, lower_gray, upper_gray) # Threshold the HSV image to get only gray colors
    ret, thresh = cv.threshold(imgray, 127, 255, cv.THRESH_BINARY_INV)
    contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_NONE)
    contours_info = [] 
    for i, c in enumerate(contours):
        max_point = c.reshape(c.shape[0], 2)[:,1].max()
        contours_info.append({"size": c.shape[0], "max_point": max_point, "contour": c})
    contours_info.sort(key=lambda x: x['size'], reverse=True)
    traces_contour = contours_info[:subfigures_num]
    traces_contour.sort(key= lambda x: x['max_point'])

    # Create traces
    for i, trace_contour in enumerate(traces_contour):
        component_name = figure_info["subfigures"][i]["text"]
        print(f"Processing trace {component_name} ...")
        trace_df = __extract_trace(
            im.shape[0],
            trace_contour["contour"],
            'a',
            component_name,
            axis_contour["contour"],
            axis_y_ticks,
            ticks_down,
            prominence_ticks,
            ignore_first_tick,
            im.copy()
        )
        traces.append(trace_df)
    return traces
