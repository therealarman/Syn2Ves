import matplotlib.pyplot as plt
import pyvista as pv
from pyvista import examples
from pyvista.core.pointset import PolyData
from pyvista.utilities.reader import STLReader
from typing import List, Tuple
import numpy as np
from PIL import Image
import cv2
import math
import time
from surface_area import rotateObj, rotate_get_scale, getWhiteVals, maxmin_tuple, itr_rotate, find_scale
from align_projection import rotate, angle, bm1


y_min, y_max, y_step = -90, 135, 45
z_min, z_max, z_step = -90, 135, 45

synapse_paths = ["./paired synapse-vesicle meshes/1 - Synapse 408.stl",
                "./paired synapse-vesicle meshes/2 - Synapse 492.stl",
                # "./paired synapse-vesicle meshes/3 - Synapse 3.stl",
                "./paired synapse-vesicle meshes/4 - Synapse 346.stl",
                "./paired synapse-vesicle meshes/5 - Synapse 473.stl",
                "./paired synapse-vesicle meshes/6 - Synapse 471.stl",
                "./paired synapse-vesicle meshes/7 - Synapse 206.stl",
                "./paired synapse-vesicle meshes/8 - Synapse 552.stl",
                "./paired synapse-vesicle meshes/9 - Synapse 167.stl",
                "./paired synapse-vesicle meshes/10 - Synapse 96.stl"]

# vesicle_path = "./paired synapse-vesicle meshes/1- Vesicle 631.stl"

'''
FIRST LOOP ITERATION

Moves virtual camera around Synapse in a circle to find angle with most viewable surface areau
'''
for synapse_path in synapse_paths:
    start_time = time.time()
    camScale = find_scale(y_min, y_max, y_step, z_min, z_max, z_step, synapse_path)
    image_array, minMaxVal = itr_rotate(y_min, y_max, y_step, z_min, z_max, z_step, synapse_path, camScale)
    print("--- %s seconds ---" % (time.time() - start_time))

    # degree_nums = range(-90, 135, 45)

    # for idx, x in enumerate(image_array):
    #     whites = getWhiteVals(x)
    #     plt.plot(degree_nums, whites, label=str(range(-90, 135, 45)[idx]))
    #     _max, _min = maxmin_tuple(degree_nums, whites)
    #     plt.scatter(_max[0], _max[1])

    # plt.legend(loc="upper right")
    # plt.xticks(np.arange(-90, 135, 45))
    # plt.show()

'''
SECOND LOOP ITERATION

This loop is fine tuned to the results of the first loop, with more steps

y_deg = minMaxVal[0]
z_deg = minMaxVal[1]

y_min, y_max = y_deg - 45, y_deg + 45
z_min, z_max = z_deg - 45, z_deg + 45

y_step = 10
z_step = 10

start_time = time.time()
hq_image_array, hq_minMaxVal = itr_rotate(y_min, y_max, y_step, z_min, z_max, z_step, synapse_path, camScale)
print("--- %s seconds ---" % (time.time() - start_time))

degree_nums = range(z_deg - 45, z_deg + 45, 10)

for idx, x in enumerate(hq_image_array):
    whites = getWhiteVals(x)
    plt.plot(degree_nums, whites, label=str(range(y_min, y_max, y_step)[idx]))
    _max, _min = maxmin_tuple(degree_nums, whites)
    plt.scatter(_max[0], _max[1])

plt.legend(loc="upper right")
plt.xticks(np.arange(z_deg - 45, z_deg + 45, 10))

plt.show()
'''