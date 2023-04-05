import matplotlib.pyplot as plt
import pyvista as pv
from typing import List, Tuple
import pandas as pd
import numpy as np
from PIL import Image, ImageEnhance, ImageOps
import os
import cv2
import math
import time
from datetime import datetime
from surface_area import getWhiteVals, maxmin_tuple
import re

def destroyPlot(plotter: pv.Plotter):
    plotter.close()

    # pv.plotting._ALL_PLOTTERS.pop(plotter._id_name)
    # plotter.renderers = None
    # delattr(plotter, 'renderers')
    plotter._theme = None
    del plotter

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
    brightArr = np.empty_like(images_arr)

    width = end - start
    norm_values = (values_arr-np.min(values_arr))/(np.max(values_arr)-np.min(values_arr)) * width + start

    for i in range(len(images_arr)):
        for j in range(len(images_arr[i])):
            img = images_arr[i][j]
            bright_img = Image.fromarray(img.astype('uint8'), 'RGB')
            # bright_img = bright_img.convert('L')
            enhancer = ImageEnhance.Brightness(bright_img)
            bright_img = enhancer.enhance(norm_values[i][j])

            if(norm_values[i][j] == np.amax(norm_values)):
                bright_img = bright_img.convert('L')
                bright_img = ImageOps.colorize(bright_img, black="black", white="green")

            brightArr[i][j] = np.array(bright_img)

    brightTile = concat_tile(brightArr)

    # plt.imshow(brightTile)
    # plt.show()
    return(brightTile)

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

def surfaceAreaAngle(path: str, camScales: list, x_range: list, y_range: list, z_og: float = 0, y_og: float = 0, center_og: tuple = (0, 0, 0)):
    surfAreaImgs = []
    surfAreaVals = []

    p = pv.Plotter(off_screen=True)

    for i in x_range:

        iterImgs = []
        iterVals = []

        for j in y_range:
            p.clear()

            sy_reader = pv.get_reader(path)
            synapse = sy_reader.read()

            synapse.rotate_z(z_og, center_og, inplace=True)
            synapse.rotate_y(y_og, center_og, inplace=True)

            synapse.rotate_x(i, synapse.center, inplace=True)
            synapse.rotate_y(j, synapse.center, inplace=True)

            p.add_mesh(synapse, lighting=False)

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

    destroyPlot(p)
    
    return(surfAreaVals, surfAreaImgs)

synDir = "C:/Users/AlexisA/Documents/2023-02-22_12-54-06 - Normal/Final Synapses"
vesDir = "C:/Users/AlexisA/Documents/2023-02-22_12-54-06 - Normal/Final Vesicles"

# pairCsv = "C:/Users/AlexisA/Documents/testPairs.xlsx"
pairCsv = "C:/Users/AlexisA/Documents/synVesPairs (Autosaved).xlsx"

def Syn2Ves(synFiles: str, vesFiles: str, pairing: str):
    camVesPos = []
    sfaVesPos = []
    vesAngle = []
    intersectVals = []
    iouVals = []

    df = pd.read_excel(pairing, engine='openpyxl')
    df1 = df[['synLabel', 'vesLabel']]
    
    synFiles = num_sort(os.listdir(synFiles))
    vesFiles = num_sort(os.listdir(vesFiles))
    
    for idx, syn in enumerate(synFiles):

        vesIdx = int(df1.iloc[idx]['vesLabel'])-1
        synIdx = int(df1.iloc[idx]['synLabel'])-1

        syn = os.path.join(synDir, synFiles[synIdx])
        ves = os.path.join(vesDir, vesFiles[vesIdx])

        print(f"{synIdx+1}, {vesIdx+1}")

        sy_reader = pv.get_reader(syn)
        ves_reader = pv.get_reader(ves)

        synapse = sy_reader.read()
        vesicle = ves_reader.read()

        p = pv.Plotter(off_screen=True)

        p1 = tuple(synapse.center)
        p2 = tuple(vesicle.center)

        # ===========================================================

        center = ((p1[0]+p2[0])/2, (p1[1]+p2[1])/2, (p1[2]+p2[2])/2)
        radius = math.sqrt(pow((p1[0] - center[0]), 2) + pow((p1[1] - center[1]), 2))
        p3 = (center[0] - radius, center[1])

        opp, adj, hyp, p3, p4 = getTriangle(p1, center, p3)

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
        radius = math.sqrt(pow((p1[0] - center[0]), 2) + pow((p1[2] - center[2]), 2))
        p3 = (center[0], center[2] - radius)

        opp, adj, hyp, p3, p4 = getTriangle((p1[0], p1[2]), (center[0], center[2]), p3)

        if(p1[0] < p3[0] and p1[2] < p3[1]):
            y_arccos = arccos(adj, hyp, opp, 2)
        elif(p1[0] > p3[0] and p1[2] > p3[1]):
            y_arccos = arccos(adj, hyp, opp, 2)
        else:
            y_arccos = -arccos(adj, hyp, opp, 2)

        synRot = synRot.rotate_y(y_arccos, center, inplace=False)
        vesRot = vesRot.rotate_y(y_arccos, center, inplace=False)

        # ===========================================================

        p.add_mesh(synRot, color="Blue", opacity = 0.5, lighting=False)
        p.add_mesh(vesRot, color="Red", opacity = 0.5, lighting=False)

        p.camera.SetParallelProjection(True)
        p.camera_position = 'xy'

        synVesCamScale = p.camera.parallel_scale
        startRot = (0, y_arccos, z_arccos)

        v0 = vesRot.center
        synapse_origin = synRot.center

        destroyPlot(p)

        # ===========================================================

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

                synapse.rotate_z(z_arccos, center, inplace=True)
                synapse.rotate_y(y_arccos, center, inplace=True)

                synapse.rotate_x(i, center, inplace=True)
                synapse.rotate_y(j, center, inplace=True)

                p.add_mesh(synapse, lighting=False)

                p.camera.SetParallelProjection(True)
                p.camera_position = 'xy'

                camScales.append(p.camera.parallel_scale)

                destroyPlot(p)

        # print(max(camScales))

        start_time = time.time()
        surfAreaVals, surfAreaImgs = surfaceAreaAngle(syn, camScales, x_range, y_range, z_arccos, y_arccos, center)
        surfAreaVals = np.array(surfAreaVals)
        maxSurfRot = np.max(surfAreaVals)
        rot_x_idx = int(np.where(surfAreaVals == maxSurfRot)[0][0])
        rot_y_idx = int(np.where(surfAreaVals == maxSurfRot)[1][0])
        surface_area_x = x_range[rot_x_idx]
        surface_area_y = y_range[rot_y_idx]
        print("Iter 1: %s seconds..." % (round(time.time() - start_time, 2)))

        x_min, x_max, x_step = surface_area_x - 45, surface_area_x + 45, 15
        x_range = range(x_min, x_max, x_step)
        y_min, y_max, y_step = surface_area_y - 45, surface_area_y + 45, 15
        y_range = range(y_min, y_max, y_step)

        start_time = time.time()
        _surfAreaVals, _surfAreaImgs = surfaceAreaAngle(syn, camScales, x_range, y_range, z_arccos, y_arccos, center)
        _surfAreaVals = np.array(_surfAreaVals)
        _maxSurfRot = np.max(_surfAreaVals)
        _rot_x_idx = int(np.where(_surfAreaVals == _maxSurfRot)[0][0])
        _rot_y_idx = int(np.where(_surfAreaVals == _maxSurfRot)[1][0])
        _surface_area_x = x_range[_rot_x_idx]
        _surface_area_y = y_range[_rot_y_idx]
        print("Iter 2: %s seconds..." % (round(time.time() - start_time, 2)))

        print(_surface_area_x)
        print(_surface_area_y)

        # rotArrImg = concat_tile(_surfAreaImgs)
        # plt.imshow(rotArrImg)
        # plt.show()

        iter1Graph = brightnessGraph(surfAreaVals, surfAreaImgs)
        iter2Graph = brightnessGraph(_surfAreaVals, _surfAreaImgs)

        cv2.imwrite(f"output/BrightGraphs/ITER1_{synIdx}_{vesIdx}.png", iter1Graph)
        cv2.imwrite(f"output/BrightGraphs/ITER2_{synIdx}_{vesIdx}.png", iter2Graph)

        ves_reader = pv.get_reader(ves)
        vesicle = ves_reader.read()

        sy_reader = pv.get_reader(syn)
        synapse = sy_reader.read()

        vesicle.rotate_z(z_arccos, center, inplace=True)
        vesicle.rotate_y(y_arccos, center, inplace=True)
        vesicle.rotate_x(_surface_area_x, synapse_origin, inplace=True)
        vesicle.rotate_y(_surface_area_y, synapse_origin, inplace=True)

        synapse.rotate_z(z_arccos, center, inplace=True)
        synapse.rotate_y(y_arccos, center, inplace=True)
        synapse.rotate_x(_surface_area_x, synapse_origin, inplace=True)
        synapse.rotate_y(_surface_area_y, synapse_origin, inplace=True)

        vesPlot = pv.Plotter(off_screen=True)
        synPlot = pv.Plotter(off_screen=True)

        synPlot.add_mesh(synapse, show_scalar_bar=False, lighting=False)
        vesPlot.add_mesh(vesicle, show_scalar_bar=False, lighting=False)

        # == Flattening ==

        synPlot.camera_position = 'xy'
        synPlot.camera.SetParallelProjection(True)
        synPlot.camera.focal_point = center
        synPlot.camera.parallel_scale = synVesCamScale
        synImg = synPlot.screenshot()

        vesPlot.camera_position = 'xy'
        vesPlot.camera.SetParallelProjection(True)
        vesPlot.camera.focal_point = center
        vesPlot.camera.parallel_scale = synVesCamScale
        vesImg = vesPlot.screenshot()

        destroyPlot(synPlot)
        destroyPlot(vesPlot)

        synapseOverlay = cv2.cvtColor(synImg, cv2.COLOR_BGR2GRAY)
        vesicleOverlay = cv2.cvtColor(vesImg, cv2.COLOR_BGR2GRAY)

        # == Intersection ==

        thresh=76
        syn_bw = cv2.threshold(synapseOverlay, thresh, 255, cv2.THRESH_BINARY)[1]
        ves_bw = cv2.threshold(vesicleOverlay, thresh, 255, cv2.THRESH_BINARY)[1]

        mergedOverlay = cv2.addWeighted(syn_bw, 0.5, ves_bw, 0.5, 0)
        intersectionImg = cv2.threshold(mergedOverlay, 128, 255, cv2.THRESH_BINARY)[1]

        # allious = cv2.hconcat([syn_bw, ves_bw, mergedOverlay, intersectionImg])
        fig, ax = plt.subplots(2,2)

        ax[0,0].imshow(syn_bw)
        ax[0,0].set_title("Synapse")
        ax[0,0].axis('off')

        ax[0,1].imshow(ves_bw)
        ax[0,1].set_title("Vesicle")
        ax[0,1].axis('off')

        ax[1,0].imshow(mergedOverlay)
        ax[1,0].set_title("Union")
        ax[1,0].axis('off')

        ax[1,1].imshow(intersectionImg)
        ax[1,1].set_title("Intersection")
        ax[1,1].axis('off')

        plt.suptitle(f"Synapse {synIdx}, Vesicle {vesIdx}")
        # plt.show()
        plt.savefig(f"output/Figs/IOU_{synIdx}_{vesIdx}.png")
        plt.close()

        mask1_area = np.count_nonzero( syn_bw )
        mask2_area = np.count_nonzero( ves_bw )
        intersection = np.count_nonzero( np.logical_and( syn_bw, ves_bw ) )
        iou = intersection/(mask1_area+mask2_area-intersection)

        print(f'IOU: {iou}')

        v1 = vesicle.center
        vectorAngle = math.degrees(np.arccos((v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]) / ((v0[0]**2 + v0[1]**2 + v0[2]**2)**0.5 * (v1[0]**2 + v1[1]**2 + v1[2]**2)**0.5)))

        print(f'Vector Angle: {vectorAngle}')

        camVesPos.append(v0)
        sfaVesPos.append(v1)
        vesAngle.append(vectorAngle)
        intersectVals.append(intersection)
        iouVals.append(iou)

    # d = {'synLabel': df1(synLabel), 'vesLabel': 'CameraPosition': camVesPos, 'VesiclePosition': sfaVesPos, 'VectorAngle': vesAngle, 'Intersect': intersectVals, 'IOU': iouVals}
    df = pd.DataFrame(data=d)

    df.to_csv("output/_FullSD3_VectorRotationData.csv")


now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start Time =", current_time)

Syn2Ves(synDir, vesDir, pairCsv)

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("End Time =", current_time)