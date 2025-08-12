# -*- coding: utf-8 -*-
"""
Modified 05/2024

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   Delete_management,
                   SpatialReference,
                   Point,
                   PointGeometry,
                   DeleteRows_management,
                   Project_management,
                   SetParameterAsText)
from arcpy.sa import (CostPath,
                      Times,
                      StreamToFeature)
from arcpy.da import (SearchCursor as SearchCursorda,
                      InsertCursor as InsertCursorda)
import os

project_name = GetParameterAsText(0)
longitude = GetParameterAsText(1)
latitude = GetParameterAsText(2)
clear_flowpaths = GetParameterAsText(3) == 'true'

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

project_folder = os.path.join(directory, "projects", project_name)

env.overwriteOutput = True
env.scratchWorkspace = project_folder
env.workspace = project_folder
env.snapRaster = os.path.join(project_folder, "dem_clip.tif")
env.extent = os.path.join(project_folder, "dem_clip.tif")

def transform_coordinates(x1, y1, sr1, sr2):
    pt1 = Point(x1, y1)
    ptgeo1 = PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    x2 = float(pt2.X)
    y2 = float(pt2.Y)
    return x2, y2

flowpath_exists = True
flowpath_index = 1
while flowpath_exists:
    line_shp = os.path.join(project_folder, "flowpath" + str(flowpath_index) +".shp")
    flowpath_exists = os.path.isfile(line_shp)
    if flowpath_exists:
        flowpath_index += 1

if clear_flowpaths:
    DeleteRows_management(os.path.join(project_folder, "addasstreams.shp"))
    for index in range(flowpath_index):
        try:
            Delete_management(os.path.join(project_folder, "flowpath" + str(index+1) +".shp"))
            Delete_management(os.path.join(project_folder, "flowpath" + str(index+1) +"_proj.shp"))
        except Exception as e:
            print(f"An error occurred: {e}")
    create_flowpath_flag = False
    output_layer = os.path.join(project_folder, "addasstreams.shp")

else:

    sr_map = SpatialReference(4326)
    sr_md = SpatialReference(26985)
    xmd, ymd = transform_coordinates(float(longitude), float(latitude), sr_map, sr_md)

    point = Point(xmd, ymd)
    ptGeometry = PointGeometry(point)
    theLine = CostPath(ptGeometry, os.path.join(project_folder, "dem_clip.tif"), os.path.join(project_folder, "flowdir.tif"), "BEST_SINGLE")
    theLine_masked = Times(theLine, os.path.join(project_folder, "basingrid.tif"))
    StreamToFeature(theLine_masked, os.path.join(project_folder, "flowdir.tif"), os.path.join(project_folder, "flowpath" + str(flowpath_index) +".shp"), "NO_SIMPLIFY")

    search_cursor = SearchCursorda(os.path.join(project_folder, "flowpath" + str(flowpath_index) +".shp"), ("OID@", "SHAPE@"))
    first_row = search_cursor.next()
    insert_cursor = InsertCursorda(os.path.join(project_folder, "addasstreams.shp"), ("OID@", "SHAPE@"))
    insert_cursor.insertRow(first_row)

    Project_management(os.path.join(project_folder, "flowpath" + str(flowpath_index) + ".shp"), os.path.join(project_folder, "flowpath" + str(flowpath_index) + "_proj.shp"), sr_map)
    create_flowpath_flag = True
    output_layer = os.path.join(project_folder, "flowpath" + str(flowpath_index) + "_proj.shp")
    
SetParameterAsText(4, create_flowpath_flag)
SetParameterAsText(5, output_layer)
