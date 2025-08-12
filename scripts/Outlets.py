# -*- coding: utf-8 -*-
"""
Modified 02/2021

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   DeleteRows_management,
                   SpatialReference,
                   Point,
                   PointGeometry,
                   RasterToPoint_conversion,
                   Delete_management,
                   management,
                   SetParameterAsText)
from arcpy.sa import (SnapPourPoint,
                      Raster)
from arcpy.da import InsertCursor
from arcpy.da import SearchCursor as SearchCursorda
import os

project_name = GetParameterAsText(0)
longitude = float(GetParameterAsText(1))
latitude = float(GetParameterAsText(2))
clear_outlets = GetParameterAsText(3) == 'true'

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

opt_folder = os.path.join(directory, "projects", project_name)

env.overwriteOutput = True
env.scratchWorkspace = opt_folder
env.workspace = opt_folder
env.snapRaster = os.path.join(opt_folder, "dem_clip.tif")

sr_map = SpatialReference(4326)
sr_md = SpatialReference(26985)

def transform_coordinates(x1, y1, sr1, sr2):
    pt1 = Point(x1, y1)
    ptgeo1 = PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    x2 = float(pt2.X)
    y2 = float(pt2.Y)
    return x2, y2

clear_outlets_flag = False
if clear_outlets:
    DeleteRows_management(os.path.join(opt_folder, "addasoutlets.shp"))

else:
    xmd, ymd = transform_coordinates(longitude, latitude, sr_map, sr_md)
    raster = Raster(os.path.join(opt_folder, "infstreams.tif"))
    cell_size = raster.meanCellWidth
    point = PointGeometry(Point(xmd, ymd), sr_md)
    outlet_snap = SnapPourPoint(point, raster, cell_size)
    outlet_point = RasterToPoint_conversion(outlet_snap, os.path.join(opt_folder, "outletSnap.shp"))
    for row in SearchCursorda(outlet_point, ["SHAPE@XY"]):
        x_new, y_new = row[0]
    xy_coordinates = (x_new, y_new)
    outlet_xy = management.GetCellValue(os.path.join(opt_folder, "infstreams.tif"), "{} {}".format(x_new, y_new))
    if outlet_xy.getOutput(0).isnumeric():
        cursor = InsertCursor(os.path.join(opt_folder, "addasoutlets.shp"), ("SHAPE@XY"))
        cursor.insertRow([xy_coordinates])
        clear_outlets_flag = True

    try:
        Delete_management(os.path.join(opt_folder, "outletxy.tif"), "")
        Delete_management(os.path.join(opt_folder, "outletSnap.shp"),"")
    except Exception as e:
        print(f"An error occurred: {e}")

SetParameterAsText(4, clear_outlets_flag)