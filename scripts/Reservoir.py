# -*- coding: utf-8 -*-
"""
Modified 02/2021

@author: Javier.Mardones
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:49:54 2020

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   SpatialReference,
                   Point,
                   PointGeometry,
                   DeleteFeatures_management,
                   GetCellValue_management,
                   SearchCursor,
                   Intersect_analysis,
                   Delete_management,
                   SpatialJoin_analysis,
                   SetParameterAsText)
from arcpy.da import InsertCursor as InsertCursorda
import os

set_i = 0
projectname = GetParameterAsText(set_i); set_i = set_i + 1
getlon = GetParameterAsText(set_i); set_i = set_i + 1
getlat = GetParameterAsText(set_i); set_i = set_i + 1

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

optfolder = os.path.join(directory, "projects", projectname)

# Environmental variables

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

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


### ADD RESERVOIR (STRUCTURE) FOR WINTR20 INPUT
xmd, ymd = coord_transf(getlon, getlat, sr_map, sr_md)
xy = (xmd, ymd)
DeleteFeatures_management(os.path.join(optfolder, "AddasReservoir.shp"))
cursor = InsertCursorda(os.path.join(optfolder, "AddasReservoir.shp"), ("SHAPE@XY"))
cursor.insertRow([xy])

in_target = os.path.join(optfolder, "AddasReservoir.shp")
in_join = os.path.join(optfolder, "subshed.shp")
out_feature_class = os.path.join(optfolder, "intersect_aux.shp")
SpatialJoin_analysis(in_target, in_join, out_feature_class)
sr = SearchCursor(out_feature_class, "", "", "ARCID", "")
for node in sr:
    rating = int(node.getValue("ARCID"))

### TRAP ERROR IN CASE TRANSECT IS IN INCORRECT REACH
arcidnode = []
fromnode = []
tonode = []
subriver_prop = SearchCursor(os.path.join(optfolder, "subrivers.shp"), "", "", "ARCID;From_Node;To_Node", "")
for node in subriver_prop:
    arcidnode.append(int(node.getValue("ARCID")))
    fromnode.append(int(node.getValue("From_Node")))
    tonode.append(int(node.getValue("To_Node")))
subreach_lst = [fn for fn in fromnode if fn in tonode]
arcid_list = []
for sub in subreach_lst:
    index = fromnode.index(sub)
    arcid_list.append(arcidnode[index])
res_trun = [os.path.join(optfolder, "AddasReservoir.shp"), os.path.join(optfolder, "subshed.shp")]
Intersect_analysis(res_trun, os.path.join(optfolder, "res_trun.shp"), "ALL", "#", "INPUT")
transect_prop = SearchCursor(os.path.join(optfolder, "res_trun.shp"), "", "", "ARCID", "")
transnodes = []
for tra in transect_prop:
    transnodes.append(int(tra.getValue("ARCID")))
condition_in = [tn for tn in transnodes if tn in arcid_list]

if len(condition_in) == 0:
    SetParameterAsText(set_i, False); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1

else:
    reselev = GetCellValue_management(os.path.join(optfolder, "dem_clip.tif"), '%s %s' % (xmd,ymd))
    reservoir_elev = str(round(float(reselev.getOutput(0)), 2))

    #Default name for reservoir
    ratingtype = "Struct"

    Delete_management(os.path.join(optfolder, "intersect_aux.shp"))
    Delete_management(os.path.join(optfolder, "xline_reach" + str(rating) + ".shp"))
    Delete_management(os.path.join(optfolder, "xpoint_reach" + str(rating) + ".shp"))

    SetParameterAsText(set_i, True); set_i = set_i + 1
    SetParameterAsText(set_i, reservoir_elev); set_i = set_i + 1
    SetParameterAsText(set_i, ratingtype); set_i = set_i + 1
    SetParameterAsText(set_i, rating); set_i = set_i + 1

Delete_management(os.path.join(optfolder, "res_trun.shp"), "")
Delete_management(os.path.join(optfolder, "intersect_aux.shp"), "")

