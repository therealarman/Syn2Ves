import matplotlib.pyplot as plt
import pyvista as pv
from pyvista import examples
from pyvista.core.pointset import PolyData
from pyvista.utilities.reader import STLReader
from typing import List, Tuple
import pandas as pd
import numpy as np
from PIL import Image, ImageEnhance, ImageDraw, ImageFont
import os
import cv2
import math
import time
import copy
from surface_area import rotateObj, rotate_get_scale, getWhiteVals, maxmin_tuple, itr_rotate, find_scale
from align_projection import rotate, angle, bm1
import re

def num_sort(sortList):
    sortList.sort(key=lambda test_string : list(map(int, re.findall(r'\d+', test_string)))[0])
    return(sortList)

def brightnessLines(images_arr, low, high, step, _low = [], _high = [], _step = []):
    degree_nums = range(low, high, step)

    if _low == []:
        _low = low
    if _high == []:
        _high = high
    if _step == []:
        _step = step

    for idx, x in enumerate(images_arr):
        whites = getWhiteVals(x)
        plt.plot(degree_nums, whites, label=str(range(_low, _high, _step)[idx]))
        _max, _min = maxmin_tuple(degree_nums, whites)
        plt.scatter(_max[0], _max[1])
    plt.legend(loc="upper right")
    plt.xticks(np.arange(-90, 135, 45))
    plt.show()

def brightnessGraph(values_arr, images_arr, start = 0.25, end = 1):
    # values_arr = np.reshape(values_arr, (-1, len(images_arr[0])))
    brightArr = np.empty_like(images_arr)

    width = end - start
    norm_values = (values_arr-np.min(values_arr))/(np.max(values_arr)-np.min(values_arr)) * width + start

    for i in range(len(images_arr)):
        for j in range(len(images_arr[i])):
            img = images_arr[i][j]
            bright_img = Image.fromarray(img.astype('uint8'), 'RGB')
            enhancer = ImageEnhance.Brightness(bright_img)
            bright_img = enhancer.enhance(norm_values[i][j])
            brightArr[i][j] = np.array(bright_img)
    brightTile = concat_tile(brightArr)

    plt.imshow(brightTile)
    plt.show()

def concat_tile(im_list_2d):
    return cv2.vconcat([cv2.hconcat(im_list_h) for im_list_h in im_list_2d])

def arccos(adj, hyp, opp, mult = 1):
    return(mult * 57.2958 * np.arccos(((adj**2 + hyp**2) - opp**2) / (2 * adj * hyp)))

def getTriangle(p1, center, p3): 
    p4 = ((p1[0]+p3[0])/2, (p1[1]+p3[1])/2)

    opp = math.dist((p1[0], p1[1]), p4)
    adj = math.dist((center[0], center[1]), p4)
    hyp = math.dist((center[0], center[1]), (p1[0], p1[1]))

    return(opp, adj, hyp, p3, p4)

def surfaceAreaAngle(path: str, x_range: list, y_range: list, z_og: float = 0, y_og: float = 0, center_og: tuple = (0, 0, 0)):
    surfAreaImgs = []
    surfAreaVals = []

    for i in x_range:

        iterImgs = []
        iterVals = []

        for j in y_range:
            sy_reader = pv.get_reader(path)
            synapse = sy_reader.read()
            p = pv.Plotter(off_screen=True)

            synapse.rotate_z(z_og, center_og, inplace=True)
            synapse.rotate_y(y_og, center_og, inplace=True)

            synapse.rotate_x(i, synapse.center, inplace=True)
            synapse.rotate_y(j, synapse.center, inplace=True)

            p.add_mesh(synapse)

            p.camera.SetParallelProjection(True)
            p.camera_position = 'xy'
            p.camera.parallel_scale = max(camScales)

            ss = p.screenshot()

            iterImgs.append(ss)

            _ss = cv2.cvtColor(ss, cv2.COLOR_BGR2GRAY)
            ret, binary = cv2.threshold(_ss,76,255,cv2.THRESH_BINARY)
            pixel_count = _ss.shape[0] * _ss.shape[1]
            white_count = np.sum(binary == 255)
            whiteRatio = round(100*white_count / pixel_count, 2)
            iterVals.append(whiteRatio)
        
        surfAreaImgs.append(iterImgs)
        surfAreaVals.append(iterVals)
    
    return(surfAreaVals, surfAreaImgs)

camVesPos = []
sfaVesPos = []
vesAngle = []

synDir = "paired synapse-vesicle meshes/Test Pairs/Synapse"
vesDir = "paired synapse-vesicle meshes/Test Pairs/Vesicle"

synFiles = num_sort(os.listdir(synDir))
vesFiles = num_sort(os.listdir(vesDir))

for idx, syn in enumerate(synFiles):
    syn = os.path.join(synDir, syn)
    ves = os.path.join(vesDir, vesFiles[idx])

    print(f"{syn}, {ves}")

    sy_reader = pv.get_reader(syn)
    ves_reader = pv.get_reader(ves)

    synapse = sy_reader.read()
    vesicle = ves_reader.read()

    p = pv.Plotter(off_screen=False)

    p1 = tuple(synapse.center)
    p2 = tuple(vesicle.center)

    center = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2)

    x_min, x_max, x_step = -90, 135, 45
    x_range = range(x_min, x_max, x_step)
    y_min, y_max, y_step = -90, 135, 45
    y_range = range(y_min, y_max, y_step)
    
    camScales = []

    for i in x_range:
        for j in y_range:
            sy_reader = pv.get_reader(syn)
            synapse = sy_reader.read()
            p = pv.Plotter(off_screen=True)

            synapse.rotate_x(i, center, inplace=True)
            synapse.rotate_y(j, center, inplace=True)

            p.add_mesh(synapse)

            p.camera.SetParallelProjection(True)
            p.camera_position = 'xy'

            camScales.append(p.camera.parallel_scale)

    print(max(camScales))

    start_time = time.time()
    surfAreaVals, surfAreaImgs = surfaceAreaAngle(syn, x_range, y_range, z_arccos, y_arccos, center)
    surfAreaVals = np.array(surfAreaVals)
    maxSurfRot = np.max(surfAreaVals)
    rot_x_idx = int(np.where(surfAreaVals == maxSurfRot)[0][0])
    rot_y_idx = int(np.where(surfAreaVals == maxSurfRot)[1][0])
    surface_area_x = x_range[rot_x_idx]
    surface_area_y = y_range[rot_y_idx]
    print("--- %s seconds ---" % (time.time() - start_time))

    x_min, x_max, x_step = surface_area_x - 45, surface_area_x + 45, 15
    x_range = range(x_min, x_max, x_step)
    y_min, y_max, y_step = surface_area_y - 45, surface_area_y + 45, 15
    y_range = range(y_min, y_max, y_step)

    start_time = time.time()
    _surfAreaVals, _surfAreaImgs = surfaceAreaAngle(syn, x_range, y_range)
    _surfAreaVals = np.array(_surfAreaVals)
    _maxSurfRot = np.max(_surfAreaVals)
    _rot_x_idx = int(np.where(_surfAreaVals == _maxSurfRot)[0][0])
    _rot_y_idx = int(np.where(_surfAreaVals == _maxSurfRot)[1][0])
    _surface_area_x = x_range[_rot_x_idx]
    _surface_area_y = y_range[_rot_y_idx]
    print("--- %s seconds ---" % (time.time() - start_time))

    print(_surface_area_x)
    print(_surface_area_y)

    ves_reader = pv.get_reader(ves)
    vesicle = ves_reader.read()
    p = pv.Plotter(off_screen=True)

    vesicle.rotate_z(z_arccos, center, inplace=True)
    vesicle.rotate_y(y_arccos, center, inplace=True)

    vesicle.rotate_x(_surface_area_x, synapse_origin, inplace=True)
    vesicle.rotate_y(_surface_area_y, synapse_origin, inplace=True)

    v1 = vesicle.center

    print(v0)
    print(v1)

    vectorAngle = np.arccos((v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]) / ((v0[0]**2 + v0[1]**2 + v0[2]**2)**0.5 * (v1[0]**2 + v1[1]**2 + v1[2]**2)**0.5))

    print(vectorAngle)

    camVesPos.append(v0)
    sfaVesPos.append(v1)
    vesAngle.append(vectorAngle)

d = {'CameraPosition': camVesPos, 'VesiclePosition': sfaVesPos, 'VectorAngle': vesAngle}
df = pd.DataFrame(data=d)

df.to_csv("output/VectorRotationData.csv")