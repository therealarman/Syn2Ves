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
    values_arr = np.reshape(values_arr, (-1, len(images_arr[0])))
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

def getTriangle(p1, center):
    radius = math.sqrt(pow((p1[0] - center[0]), 2) + pow((p1[1] - center[1]), 2))

    p3 = (center[0] - radius, center[1])
    p4 = ((p1[0]+p3[0])/2, (p1[1]+p3[1])/2)

    opp = math.dist((p1[0], p1[1]), p4)
    adj = math.dist((center[0], center[1]), p4)
    hyp = math.dist((center[0], center[1]), (p1[0], p1[1]))

    return(opp, adj, hyp, p3, p4)

y_min, y_max, y_step = -90, 135, 45
z_min, z_max, z_step = -90, 135, 45

# synapse_paths = ["./paired synapse-vesicle meshes/1 - Synapse 408.stl",
#                 # "./paired synapse-vesicle meshes/2 - Synapse 492.stl",
#                 "./paired synapse-vesicle meshes/3 - Synapse 3.stl",
#                 "./paired synapse-vesicle meshes/4 - Synapse 346.stl",
#                 "./paired synapse-vesicle meshes/5 - Synapse 473.stl",
#                 "./paired synapse-vesicle meshes/6 - Synapse 471.stl",
#                 "./paired synapse-vesicle meshes/7 - Synapse 206.stl",
#                 "./paired synapse-vesicle meshes/8 - Synapse 552.stl",
#                 "./paired synapse-vesicle meshes/9 - Synapse 167.stl",
#                 "./paired synapse-vesicle meshes/10 - Synapse 96.stl"]

# vesicle_paths = ["./paired synapse-vesicle meshes/1 - Vesicle 631.stl",
#                 # "./paired synapse-vesicle meshes/2 - Vesicle 698.stl",
#                 "./paired synapse-vesicle meshes/3 - Vesicle 284.stl",
#                 "./paired synapse-vesicle meshes/4 - Vesicle 275.stl",
#                 "./paired synapse-vesicle meshes/5 - Vesicle 480.stl",
#                 "./paired synapse-vesicle meshes/6 - Vesicle 290.stl",
#                 "./paired synapse-vesicle meshes/7 - Vesicle 604.stl",
#                 "./paired synapse-vesicle meshes/8 - Vesicle 752.stl",
#                 "./paired synapse-vesicle meshes/9 - Vesicle 182.stl",
#                 "./paired synapse-vesicle meshes/10 - Vesicle 66.stl"]

SurfaceAreaRots = []
CameraAlignRots = []
Difference = []
zDiffs = []
yDiffs = []
IntersectVals = []
IOUs = []

synDir = "paired synapse-vesicle meshes/Test Pairs/Synapse"
vesDir = "paired synapse-vesicle meshes/Test Pairs/Vesicle"

synFiles = num_sort(os.listdir(synDir))
vesFiles = num_sort(os.listdir(vesDir))

for idx, syn in enumerate([synFiles[0]]):
# for idx, syn in enumerate(synFiles):

    syn = os.path.join(synDir, syn)
    ves = os.path.join(vesDir, vesFiles[idx])

    print(f"{syn}, {ves}")
    
    # start_time = time.time()
    camScale = find_scale(y_min, y_max, y_step, z_min, z_max, z_step, syn)
    image_array, minMaxVal, whiteArr = itr_rotate(y_min, y_max, y_step, z_min, z_max, z_step, syn, camScale)
    # print("--- %s seconds ---" % (time.time() - start_time))
    # brightnessGraph(whiteArr, image_array, 0.1, 1)

    sy_reader = pv.get_reader(syn)
    ves_reader = pv.get_reader(ves)

    synapse = sy_reader.read()
    vesicle = ves_reader.read()

    p = pv.Plotter(off_screen=True)

    p1 = tuple(synapse.center)
    p2 = tuple(vesicle.center)

    # ===========================================================

    center = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2)
    opp, adj, hyp, p3, p4 = getTriangle(p1, center)

    if(p1[1] > p3[1]):
        z_arccos = arccos(adj, hyp, opp, 2)
    else:
        z_arccos = -arccos(adj, hyp, opp, 2)

    synRot = synapse.rotate_z(z_arccos, center, inplace=False)
    vesRot = vesicle.rotate_z(z_arccos, center, inplace=False)

    # ===========================================================

    p1 = tuple(synRot.center)
    p2 = tuple(vesRot.center)

    center = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2)
    opp, adj, hyp, p3, p4 = getTriangle((p1[0], p1[2]), (center[0], center[2]))

    if(p1[0] < p3[0] and p1[2] < p3[1]):
        y_arccos = -arccos(adj, hyp, opp, 2)
    elif(p1[0] > p3[0] and p1[2] > p3[1]):
        y_arccos = -arccos(adj, hyp, opp, 2)
    else:
        y_arccos = arccos(adj, hyp, opp, 2)

    synRot = synRot.rotate_y(y_arccos, center, inplace=False)
    vesRot = vesRot.rotate_y(y_arccos, center, inplace=False)

    # ===========================================================

    # p.close()

    p.add_mesh(synRot, color="Blue", opacity = 0.5)
    p.add_mesh(vesRot, color="Red", opacity = 0.5)

    p.camera.SetParallelProjection(True)
    p.camera_position = 'yz'
    camScale = p.camera.parallel_scale

    #cv2.imwrite(f"output\{idx} - CamAlign.png", p.screenshot())

    # p.close()

    p.show()

    # == Surface Area ==

    p = pv.Plotter(off_screen=True)

    synapseRot = copy.copy(synapse)
    vesicleRot = copy.copy(vesicle)

    synapseRot.rotate_y(minMaxVal[0], synapse.center, inplace=True)
    synapseRot.rotate_z(minMaxVal[1], synapse.center, inplace=True)

    vesicleRot.rotate_y(minMaxVal[0], synapse.center, inplace=True)
    vesicleRot.rotate_z(minMaxVal[1], synapse.center, inplace=True)

    normal = (1, 0, 0)
    sy_projected = synapseRot.project_points_to_plane(origin=center, normal=normal)
    ves_projected = vesicleRot.project_points_to_plane(origin=center, normal=normal)

    p.add_mesh(synapseRot, color="Blue")
    p.add_mesh(vesicleRot, color="Red")

    p.camera_position = 'yz'
    p.camera.SetParallelProjection(True)
    p.camera.focal_point = center
    camScale = p.camera.parallel_scale

    #cv2.imwrite(f"output\{idx} - SurfAlign.png", p.screenshot())

    p.close()

    synPlot = pv.Plotter(off_screen=True)
    vesPlot = pv.Plotter(off_screen=True)

    synPlot.add_mesh(sy_projected, show_scalar_bar=False)
    vesPlot.add_mesh(ves_projected, show_scalar_bar=False)

    # == Flattening ==

    synPlot.camera_position = 'yz'
    synPlot.camera.SetParallelProjection(True)
    synPlot.camera.focal_point = center
    synPlot.camera.parallel_scale = camScale
    synImg = synPlot.screenshot()

    vesPlot.camera_position = 'yz'
    vesPlot.camera.SetParallelProjection(True)
    vesPlot.camera.focal_point = center
    vesPlot.camera.parallel_scale = camScale
    vesImg = vesPlot.screenshot()

    synapseOverlay = cv2.cvtColor(synImg, cv2.COLOR_BGR2GRAY)
    vesicleOverlay = cv2.cvtColor(vesImg, cv2.COLOR_BGR2GRAY)

    # == Intersection ==

    thresh=76
    syn_bw = cv2.threshold(synapseOverlay, thresh, 255, cv2.THRESH_BINARY)[1]
    ves_bw = cv2.threshold(vesicleOverlay, thresh, 255, cv2.THRESH_BINARY)[1]

    mergedOverlay = cv2.addWeighted(syn_bw, 0.5, ves_bw, 0.5, 0)
    intersectionImg = cv2.threshold(mergedOverlay, 128, 255, cv2.THRESH_BINARY)[1]

    #cv2.imwrite(f"output\{idx} - Intersect.png", intersectionImg)

    # intersection = (100 * np.sum(intersectionImg == 255)) / (intersectionImg.shape[0] * intersectionImg.shape[1])
    # print(f"{intersectionImg.shape[0] * intersectionImg.shape[1]}, {np.sum(intersectionImg == 255)}")

    mask1_area = np.count_nonzero( syn_bw )
    mask2_area = np.count_nonzero( ves_bw )
    intersection = np.count_nonzero( np.logical_and( syn_bw, ves_bw ) )
    iou = intersection/(mask1_area+mask2_area-intersection)

    # == Vals ==

    z_diff = z_arccos - minMaxVal[1]
    y_diff = y_arccos - minMaxVal[0]

    zDiffs.append(z_diff)
    yDiffs.append(y_diff)

    diff = (abs(z_diff) + abs(y_diff)) / 2

    SurfaceAreaRots.append( (z_arccos, y_arccos) )
    CameraAlignRots.append( (minMaxVal[1], minMaxVal[0]) )
    Difference.append( diff )
    IntersectVals.append(intersection)
    IOUs.append(iou)

d = {'SurfArea': SurfaceAreaRots, 'CamAlign': CameraAlignRots, 'Z Diff': zDiffs, 'Y Diff': yDiffs, 'Diff': Difference, 'Intersect': IntersectVals, 'IOU': IOUs}
df = pd.DataFrame(data=d)

df.to_csv("output/RotationData.csv")