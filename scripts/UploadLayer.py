# -*- coding: utf-8 -*-
"""
Modified 05/2024

@author: Javier.Mardones
"""

import arcpy
import os, ast

def get_min_elevation(filepath):
    min_elev = float('inf')
    with arcpy.da.SearchCursor(filepath, "Elev") as cursor:
        for row in cursor:
            if float(row[0]) < min_elev:
                min_elev = float(row[0])
    return min_elev

def delete_non_min_elevation(filepath, min_elev):
    with arcpy.da.UpdateCursor(filepath, "Elev") as cursor:
        for row in cursor:
            if not float(row[0]) == min_elev:
                cursor.deleteRow()


projectname = arcpy.GetParameterAsText(0)
incoords = arcpy.GetParameterAsText(1)
cellsize = arcpy.GetParameterAsText(2)

# get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

optfolder = os.path.join(directory, "projects")
optfolder = os.path.join(optfolder, projectname)

arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = optfolder
arcpy.env.workspace = optfolder

sr_map = arcpy.SpatialReference(4326)
sr_md = arcpy.SpatialReference(26985)

in_json_file = os.path.join(optfolder, "uploaded_layer.txt")

with open(in_json_file, 'w') as f:
    f.write(incoords)

out_features = os.path.join(optfolder, "uploaded_layer_proj.shp")

incoords = ast.literal_eval(incoords)

array = arcpy.Array()
for x in incoords[0]:
    array.add(arcpy.Point(float(x[0]), float(x[1])))

m = arcpy.Polygon(array, sr_map)
arcpy.CopyFeatures_management(m, out_features)

arcpy.Project_management(out_features, os.path.join(
    optfolder, "uploaded_layer.shp"), sr_md)

arcpy.CopyFeatures_management(os.path.join(optfolder, "uploaded_layer.shp"), os.path.join(optfolder, "wshed.shp"))
arcpy.PolygonToRaster_conversion(os.path.join(optfolder, "wshed.shp"), "FID", os.path.join(optfolder, "basingrid.tif"), "", "", cellsize)
arcpy.GeneratePointsAlongLines_management(os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "usr_pour_point.shp"), "Distance", cellsize)
arcpy.sa.ExtractMultiValuesToPoints(os.path.join(optfolder, "usr_pour_point.shp"), [[os.path.join(optfolder, "dem_clip.tif"), "Elev"]])

min_elev = get_min_elevation(os.path.join(optfolder, "usr_pour_point.shp"))
delete_non_min_elevation(os.path.join(optfolder, "usr_pour_point.shp"), min_elev)
