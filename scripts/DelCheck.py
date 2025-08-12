# -*- coding: utf-8 -*-
"""
Modified 05/2024

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   Extent,
                   Point,
                   PointGeometry,
                   Raster,
                   RasterToPoint_conversion,
                   SpatialReference,
                   management,
                   SetParameterAsText)
from arcpy.sa import SnapPourPoint
from arcpy.da import SearchCursor as SearchCursorda
import os
import json

GP_TASK_VERSION = "0.96"

mouse_lat_proj = float(GetParameterAsText(0));
mouse_lon_proj = float(GetParameterAsText(1));
cellsize = int(GetParameterAsText(2));
userstream = GetParameterAsText(3) == "true";
projectname = GetParameterAsText(4);

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

if userstream:
    optfolder = os.path.join(directory, "projects")
    optfolder = os.path.join(optfolder, projectname)
    infstreams = os.path.join(optfolder, "infstr_ext.tif")
else:
    infstreams = os.path.join(directory, "data","hydro{}".format(cellsize) + ".gdb", "infstreams{}".format(cellsize))

try:
    versionfile = os.path.join(directory, "data", "dataversion.txt")
    with open(versionfile) as f:
        dataversion = f.readlines()
except:
    dataversion = ["-"]

sr_map = SpatialReference(4326)
sr_md = SpatialReference(26985)

def coord_transf(x1,y1,sr1,sr2):
    pt1 = Point(x1,y1)
    ptgeo1 = PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    x2 = float(pt2.X)
    y2 = float(pt2.Y)
    return x2, y2

x, y = coord_transf(mouse_lon_proj, mouse_lat_proj, sr_map, sr_md)
rast = Raster(infstreams)
cellsize = rast.meanCellWidth
point = PointGeometry(Point(x, y))

env.overwriteOutput = True
env.extent = Extent(x-cellsize, y-cellsize, x+cellsize, y+cellsize)
env.snapRaster = infstreams

outSnapPour = SnapPourPoint(point, infstreams, cellsize)

pour_point = RasterToPoint_conversion(outSnapPour)

for row in SearchCursorda(pour_point, ["SHAPE@XY"]):
    x, y = row[0]
outletxy = management.GetCellValue(infstreams, "{} {}".format(x, y))

xmap, ymap = coord_transf(x, y, sr_md, sr_map)

if outletxy.getOutput(0).isnumeric():
    SetParameterAsText(5, True)
else:
    SetParameterAsText(5, False)

SetParameterAsText(6, pour_point)
SetParameterAsText(7, json.dumps([xmap, ymap]))
SetParameterAsText(8, GP_TASK_VERSION)
SetParameterAsText(9, dataversion[0])

del pour_point, outletxy, outSnapPour




