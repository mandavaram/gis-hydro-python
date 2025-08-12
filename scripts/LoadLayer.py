# -*- coding: utf-8 -*-
"""
Modified: 05/2024

@author: Javier.Mardones
"""

import os
from arcpy import (GetParameterAsText, env, SpatialReference, RasterToPolygon_conversion, Dissolve_management,
                   Project_management, UpdateCursor, Delete_management, Merge_management, SetParameterAsText)
from arcpy.sa import Contour, Times

SR_MAP_WKID = 4326
directory = os.environ['GISHydro_DIR']
sr_map = SpatialReference(SR_MAP_WKID)

projectname = GetParameterAsText(0)
userchoice = GetParameterAsText(1)
contourInterval = float(GetParameterAsText(2))
baseContour = float(GetParameterAsText(3))
reaches = int(GetParameterAsText(4))

optfolder = os.path.join(directory, "projects", projectname)

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = os.path.join(optfolder, "dem_clip.tif")
env.extent = env.snapRaster

def contours(optfolder, sr_map, contourInterval, baseContour):
    contours = Times(os.path.join(optfolder, "dem_clip.tif"), os.path.join(optfolder, "basingrid.tif"))
    contours = Contour(contours, os.path.join(optfolder, "contours.shp"), contourInterval, baseContour)
    Project_management(contours, os.path.join(optfolder, "contours_aux.shp"), sr_map)

    cont = UpdateCursor(os.path.join(optfolder, "contours_aux.shp"))
    for c in cont:
        c.ID = 1
        cont.updateRow(c)

    Dissolve_management(os.path.join(optfolder, "contours_aux.shp"), os.path.join(optfolder, "contours_proj.shp"), "ID")
    Delete_management(os.path.join(optfolder, "contours_aux.shp"), "")
    
    return os.path.join(optfolder, "contours_proj.shp")

def longpath(optfolder, sr_map, subareas):
    paths = [os.path.join(optfolder, f"Longest_Path_Sub_{i+1}.shp") for i in range(reaches)]
    Merge_management(paths, os.path.join(optfolder, "paths_aux.shp"))
    Project_management(os.path.join(optfolder, "paths_aux.shp"), os.path.join(optfolder, "paths_proj.shp"), sr_map)

    Delete_management(os.path.join(optfolder, "paths_aux.shp"), "")
    
    return os.path.join(optfolder, "paths_proj.shp")

def landuse(optfolder, sr_map):
    RasterToPolygon_conversion(os.path.join(optfolder, "landuse.tif"), os.path.join(optfolder, "landuse_aux.shp"),"NO_SIMPLIFY","CLASS_NAME")
    Dissolve_management(os.path.join(optfolder, "landuse_aux.shp"), os.path.join(optfolder, "landuse.shp"), "CLASS_NAME")
    Project_management(os.path.join(optfolder, "landuse.shp"), os.path.join(optfolder, "landuse_proj.shp"), sr_map)

    Delete_management(os.path.join(optfolder, "landuse.shp"), "")
    Delete_management(os.path.join(optfolder, "landuse_aux.shp"), "")
    
    return os.path.join(optfolder, "landuse_proj.shp")

def soils(optfolder, sr_map):
    RasterToPolygon_conversion(os.path.join(optfolder, "soils.tif"), os.path.join(optfolder, "soils_aux.shp"),"NO_SIMPLIFY","Value")
    Dissolve_management(os.path.join(optfolder, "soils_aux.shp"), os.path.join(optfolder, "soils.shp"), "gridcode")
    Project_management(os.path.join(optfolder, "soils.shp"), os.path.join(optfolder, "soils_proj.shp"), sr_map)

    Delete_management(os.path.join(optfolder, "soils.shp"), "")
    Delete_management(os.path.join(optfolder, "soils_aux.shp"), "")
    
    return os.path.join(optfolder, "soils_proj.shp")

def curvenumber(optfolder, sr_map):
    RasterToPolygon_conversion(os.path.join(optfolder, "curvenumber.tif"), os.path.join(optfolder, "curvenumber_aux.shp"),"NO_SIMPLIFY","Value")
    Dissolve_management(os.path.join(optfolder, "curvenumber_aux.shp"), os.path.join(optfolder, "curvenumber.shp"), "gridcode")
    Project_management(os.path.join(optfolder, "curvenumber.shp"), os.path.join(optfolder, "curvenumber_proj.shp"), sr_map)

    Delete_management(os.path.join(optfolder, "curvenumber.shp"), "")
    Delete_management(os.path.join(optfolder, "curvenumber_aux.shp"), "")
    
    return os.path.join(optfolder, "curvenumber_proj.shp")

if userchoice == "Contours":
    output_layer = contours(optfolder, sr_map, contourInterval, baseContour)
elif userchoice == "Longest Path":
    output_layer = longpath(optfolder, sr_map, reaches)
elif userchoice == "Land Use":
    output_layer = landuse(optfolder, sr_map)
elif userchoice == "Soils":
    output_layer = soils(optfolder, sr_map)
elif userchoice == "Curve Number":
    output_layer = curvenumber(optfolder, sr_map)

SetParameterAsText(5, output_layer)
