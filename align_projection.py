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

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def angle(a, b, c):
    ang = math.degrees(math.atan2(c[1]-b[1], c[0]-b[0]) - math.atan2(a[1]-b[1], a[0]-b[0]))
    return ang + 360 if ang < 0 else ang

'''
def plotObj(polyObj: PolyData, focalPoint: Tuple[int, int, int], camScale: float):
    p = pv.Plotter(off_screen=True)
    
    _projected = polyObj.project_points_to_plane(origin=polyObj.center, normal=normal)
    
    p.add_mesh(_projected)
    
    p.camera_position = 'yz'
    p.camera.SetParallelProjection(True)
    p.camera.focal_point = focalPoint
    p.camera.parallel_scale = camScale
        
    p.show()
    
    return(p.screenshot())
'''

def bm1(mask1, mask2):
    mask1_area = np.count_nonzero( mask1 )
    mask2_area = np.count_nonzero( mask2 )
    intersection = np.count_nonzero( np.logical_and( mask1, mask2 ) )
    iou = intersection/(mask1_area+mask2_area-intersection)
    return iou