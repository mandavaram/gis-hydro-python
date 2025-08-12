# -*- coding: utf-8 -*-
"""
Modified 02/2021

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   SpatialReference,
                   Point,
                   PointGeometry,
                   Array,
                   Polyline,
                   CopyFeatures_management,
                   SearchCursor,
                   Intersect_analysis,
                   Delete_management,
                   CreateRandomPoints_management,
                   UpdateCursor,
                   Raster,
                   SpatialJoin_analysis,
                   SetParameterAsText)
from arcpy.sa import (ExtractMultiValuesToPoints,
                      ZonalStatisticsAsTable)
from arcpy.da import SearchCursor as SearchCursorda
import os
from operator import mul
import json

set_i = 0
projectname = GetParameterAsText(set_i); set_i = set_i + 1
getlon = json.loads(GetParameterAsText(set_i)); set_i = set_i + 1
getlat = json.loads(GetParameterAsText(set_i)); set_i = set_i + 1
nMain_val = float(GetParameterAsText(set_i)); set_i = set_i + 1
nLeft_val = float(GetParameterAsText(set_i)); set_i = set_i + 1
nRight_val = float(GetParameterAsText(set_i)); set_i = set_i + 1
user_rs = GetParameterAsText(set_i); set_i = set_i + 1
user_be = GetParameterAsText(set_i); set_i = set_i + 1
user_cw = GetParameterAsText(set_i); set_i = set_i + 1
user_cd = GetParameterAsText(set_i); set_i = set_i + 1

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])
directorygdb = os.path.join(directory, "data","gishydro.gdb")
optfolder = os.path.join(directory, "projects", projectname)

# Environmental variables

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

def coord_transf(x1,y1,sr1,sr2):
    pt1 = Point(x1,y1)
    ptgeo1 = PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    x2 = float(pt2.X)
    y2 = float(pt2.Y)
    return x2, y2

# *******************************************************************************************************
# Add DEM elevation field to transect line
# *******************************************************************************************************

sr_map = SpatialReference(4326)
sr_md = SpatialReference(26985)

point = Point()
array = Array()
for lon,lat in zip(getlon,getlat):
    xmd, ymd = coord_transf(lon, lat, sr_map, sr_md)
    point.X = xmd
    point.Y = ymd
    array.add(point)

polyline = Polyline(array)
CopyFeatures_management(polyline, os.path.join(optfolder, "trans_line.shp"))  # saving transect line

### TRAP ERROR IN CASE TRANSECT IS IN INCORRECT REACH
arcidnode = []
fromnode = []
tonode = []
subriver_prop = SearchCursor(os.path.join(optfolder, "subrivers.shp"), "", "", "ARCID;From_Node;To_Node", "")
for node in subriver_prop:
    arcidnode.append(int(node.getValue("ARCID")))
    fromnode.append(int(node.getValue("From_Node")))
    tonode.append(int(node.getValue("To_Node")))
subreach_lst = [x for x in fromnode if x in tonode]
arcid_list = []
for sub in subreach_lst:
    index = fromnode.index(sub)
    arcid_list.append(arcidnode[index])
line_trun = [os.path.join(optfolder, "trans_line.shp"), os.path.join(optfolder, "subshed.shp")]
Intersect_analysis(line_trun, os.path.join(optfolder, "line_trun.shp"), "ALL", "#", "INPUT")
transect_prop = SearchCursor(os.path.join(optfolder, "line_trun.shp"), "", "", "ARCID", "")
transnodes = []
for tra in transect_prop:
    transnodes.append(int(tra.getValue("ARCID")))
condition_in = [x for x in transnodes if x in arcid_list]

if len(condition_in) == 0:

    SetParameterAsText(set_i, False); set_i = set_i + 1 #transect not in reach
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1
    SetParameterAsText(set_i, ""); set_i = set_i + 1

    Delete_management(os.path.join(optfolder, "trans_line.shp"), "")
    Delete_management(os.path.join(optfolder, "line_trun.shp"), "")

else:

    CreateRandomPoints_management(optfolder, "transect.shp", os.path.join(optfolder, "line_trun.shp"), "#", "#", 6, "#", "#")  # 02-15-2015: changed sampling interval from 3 to 6

    # extract multivalues to points using DEM
    ExtractMultiValuesToPoints(os.path.join(optfolder, "transect.shp"), os.path.join(optfolder, "dem_clip.tif"), "NONE")
    ExtractMultiValuesToPoints(os.path.join(optfolder, "transect.shp"), os.path.join(optfolder, "flowacc.tif"), "NONE")

    # correct elevation values in feet
    rows = UpdateCursor(os.path.join(optfolder, "transect.shp"))
    for row in rows:
        row.setValue("dem_clip", row.getValue("dem_clip"))  # elevation values are already in feet -- no need for conversion
        rows.updateRow(row)

    flowacc_list = []
    upstreamDA = os.path.join(optfolder, "transect.shp")
    usda = SearchCursor(upstreamDA, "", "", "flowacc", "")
    for da in usda:
        flow = da.getValue("flowacc")
        flowacc_list.append(flow)


    # *******************************************************************************************************
    # Calculate upstream drainage area
    # *******************************************************************************************************
    upda_max = max(flowacc_list)
    flowacc_cell = os.path.join(optfolder, "flowacc.tif")
    rast = Raster(flowacc_cell)
    cellsize = rast.meanCellWidth
    cellsq = cellsize * cellsize

    areami2_usda = float((upda_max * cellsq) / 2588881)  # conversion into sq miles
    areami2_usda = "{0:.2f}".format(areami2_usda)

    # *******************************************************************************************************
    # get reach slope using subsheds, transect, and subriver shapefiles
    # *******************************************************************************************************

    if user_rs and user_rs != "#":
        reachslope = float(user_rs)
    else:
        rSlope_intersect = [os.path.join(optfolder, "line_trun.shp"), os.path.join(optfolder, "subshed.shp")]
        Intersect_analysis(rSlope_intersect, os.path.join(optfolder, "reachslope.shp"), "ALL", "#", "INPUT")
        reach = SearchCursor(os.path.join(optfolder, "reachslope.shp"), "", "", "Slope", "")
        for row in reach:
            reachslope = round(row.getValue("Slope"),5)
        del row
        del reach

    # *******************************************************************************************************
    # get bankfull channel width and depth
    # *******************************************************************************************************
    prov = os.path.join(directorygdb, "md_phyisio_provinces")
    ZonalStatisticsAsTable(prov, "PROVINCE", os.path.join(optfolder, "basingrid.tif"), "theVTab", "DATA", "ALL")

    sumarea = 0
    theVTab = os.path.join(optfolder, "theVTab.dbf")
    theVTab = SearchCursor("theVTab", "", "", "Count", "")
    for each in theVTab:
        count = each.getValue("Count")
        sumarea = sumarea + count
    sumArea = sumarea

    del each

    AParea = 0
    PDarea = 0
    CParea = 0
    theVTab = os.path.join(optfolder, "theVTab.dbf")
    theVTab = SearchCursor("theVTab", "", "", "Province;Count", "")
    for row in theVTab:
        theArea = float(row.getValue("Count"))
        areapercent = float((theArea / sumArea) * 100)
        if row.getValue("Province") == "A":
            AParea = AParea + float(areapercent)
        elif row.getValue("Province") == "B":
            AParea = AParea + float(areapercent)
        elif row.getValue("Province") == "P":
            PDarea = PDarea + float(areapercent)
        elif row.getValue("Province") == "W":
            CParea = CParea + float(areapercent)
        elif row.getValue("Province") == "E":
            CParea = CParea + float(areapercent)
    del row

    regionarea = [AParea / 100, PDarea / 100, CParea / 100]

    a = [13.87, 14.78, 10.3]
    b = [0.44, 0.39, 0.38]
    c = [0.95, 1.18, 1.01]
    d = [0.31, 0.34, 0.32]

    Coef_W = sum(map(mul, regionarea, a))
    Exp_W = sum(map(mul, regionarea, b))

    if user_cw and user_cw != "#":
        Wbf = float(user_cw)
    else:
        Wbf = round(Coef_W * (float(areami2_usda) ** Exp_W),2)

    Coef_D = sum(map(mul, regionarea, c))
    Exp_D = sum(map(mul, regionarea, d))

    if user_cd and user_cd != "#":
        Dbf = float(user_cd)
    else:
        Dbf = round(Coef_D * (float(areami2_usda) ** Exp_D),2)

                ## get elevation and flow accumulation values along transect line
    elev = []
    totlength = 0
    count = 0
    with SearchCursorda(os.path.join(optfolder, "trans_line.shp"), 'SHAPE@') as rows:
        for row in rows:
            totlength = totlength + row[0].length* 3.28084
    theelev = SearchCursor(os.path.join(optfolder, "transect.shp"), "", "", "dem_clip", "")
    for row in theelev:
        count = count + 1
        elev.append(row.getValue("dem_clip"))


    u = 0.65
    length = [x*totlength/count for x in range(count)]

    minindx = elev.index(min(elev))
    if user_be and user_be != "#":
        bfelev = float(user_be)
    else:
        bfelev = min(elev)
    alpha = (Wbf / 2) / (Dbf ** u)
    xcentr = length[minindx]
    del elev[minindx]
    del length[minindx]

    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 4)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 3)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 2)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 1)
    elev.insert(minindx, bfelev - Dbf)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 1)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 2)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 3)
    elev.insert(minindx, (bfelev - Dbf) + (Dbf / 4) * 4)

    length.insert(minindx, xcentr + (alpha * (((Dbf / 4) * 4) ** u)))
    length.insert(minindx, xcentr + (alpha * (((Dbf / 4) * 3) ** u)))
    length.insert(minindx, xcentr + (alpha * (((Dbf / 4) * 2) ** u)))
    length.insert(minindx, xcentr + (alpha * (((Dbf / 4) * 1) ** u)))
    length.insert(minindx, xcentr)
    length.insert(minindx, xcentr - (alpha * (((Dbf / 4) * 1) ** u)))
    length.insert(minindx, xcentr - (alpha * (((Dbf / 4) * 2) ** u)))
    length.insert(minindx, xcentr - (alpha * (((Dbf / 4) * 3) ** u)))
    length.insert(minindx, xcentr - (alpha * (((Dbf / 4) * 4) ** u)))

    plot_data = [length, elev]


    # ******************************************************************************************************
    # get reach slope, Wbf, Dbf, Manning numbers, reach number, stream it drains to number,
    # distance, and elevation to include in rating table input text file
    # ******************************************************************************************************

    in_target = os.path.join(optfolder, "trans_line.shp")
    in_join = os.path.join(optfolder, "subrivers.shp")
    out_feature_class = os.path.join(optfolder, "intersect_aux.shp")
    SpatialJoin_analysis(in_target, in_join, out_feature_class)
    sr = SearchCursor(out_feature_class, "", "", "ARCID", "")
    for node in sr:
        rating = int(node.getValue("ARCID"))
    
    CopyFeatures_management(polyline, os.path.join(optfolder, "xs_" + str(rating) + ".shp"))

    # ******************************************************************************************************
    # run rating table algorithm (replaced rattab.exe)
    # ******************************************************************************************************

    p  = (((Dbf / 4) * 1 - (Dbf / 4) * 0)**2 +
          (Wbf / 2 - alpha * (((Dbf / 4) * 3) ** u) - (Wbf / 2 - alpha * (((Dbf / 4) * 4) ** u)))**2)**0.5
    p += (((Dbf / 4) * 2 - (Dbf / 4) * 1)**2 +
          (Wbf / 2 - alpha * (((Dbf / 4) * 2) ** u) - (Wbf / 2 - alpha * (((Dbf / 4) * 3) ** u)))**2)**0.5
    p += (((Dbf / 4) * 3 - (Dbf / 4) * 2)**2 +
          (Wbf / 2 - alpha * (((Dbf / 4) * 1) ** u) - (Wbf / 2 - alpha * (((Dbf / 4) * 2) ** u)))**2)**0.5
    p += (((Dbf / 4) * 4 - (Dbf / 4) * 3)**2 +
          (Wbf / 2 - alpha * (((Dbf / 4) * 0) ** u) - (Wbf / 2 - alpha * (((Dbf / 4) * 1) ** u)))**2)**0.5
    pch = p*2

    n = 30  #number of iterations
    maxleft = max(elev[0:elev.index(min(elev))])
    maxright = max(elev[elev.index(min(elev))+1:])
    dx = (min(maxleft, maxright) - min(elev) + Dbf)/float(n)

    stage = [min(elev)]
    discharge = [0]
    endarea = [0]
    topwidth = [0]

    for i in range(n):
        areal = 0
        periml = 0
        arear = 0
        perimr = 0
        width = 0
        bar = min(elev) + Dbf + dx*i
        for j in range(len(elev)-1):
            xj = 0
            ajl = 0
            ajr = 0
            pjl = 0
            pjr = 0
            diff = bar - elev[j+1]
            if length[j] < length[elev.index(min(elev))]:
                if diff > 0:
                    if elev[j] >= bar:
                        yj = diff
                        xj = (length[j+1]-length[j])*yj/(elev[j]-elev[j+1])
                        ajl = yj*xj*0.5
                        pjl = (yj**2 + xj**2)**0.5
                    else:
                        yj = (diff + bar - elev[j])/2
                        xj = length[j+1] - length[j]
                        ajl = yj*xj
                        pjl = ((elev[j]-elev[j+1])**2 + (xj)**2)**0.5
                elif elev[j] < bar:
                        yj = bar - elev[j]
                        xj = (length[j+1]-length[j])*yj/(elev[j+1]-elev[j])
                        ajl = yj*xj*0.5
                        pjl = (yj**2 + xj**2)**0.5
            else:
                if diff > 0:
                    if elev[j] >= bar:
                        yj = diff
                        xj = (length[j+1]-length[j])*yj/(elev[j]-elev[j+1])
                        ajr = yj*xj*0.5
                        pjr = (yj**2 + xj**2)**0.5
                    else:
                        yj = (diff + bar - elev[j])/2
                        xj = length[j+1] - length[j]
                        ajr = yj*xj
                        pjr = ((elev[j]-elev[j+1])**2 + (xj)**2)**0.5
                elif elev[j] < bar:
                        yj = bar - elev[j]
                        xj = (length[j+1]-length[j])*yj/(elev[j+1]-elev[j])
                        ajr = yj*xj*0.5
                        pjr = (yj**2 + xj**2)**0.5
            width += xj
            areal += ajl
            periml += pjl
            arear += ajr
            perimr += pjr

        if width > 0 and bar < max(elev):
            hr = (areal + arear)/(periml + perimr)
            ncomp = ((periml * (nLeft_val)**(1.5) + pch * nMain_val**(1.5) + perimr * (nRight_val)**(1.5)) / (periml + perimr))**(0.66667) ## Equal velocity method

            stage.append(bar)
            discharge.append(1.49 * (areal + arear) * hr**(0.66667) * reachslope**0.5 / ncomp)
            endarea.append(areal + arear)
            topwidth.append(width)

            if discharge[-1] < discharge[-2]:
                discharge[-1] = discharge[-2] + 1


    discharge = [float(elem) for elem in discharge]
    discharge.sort()

    ST = ['%.1f' % elem for elem in stage]
    DI = ['%.1f' % elem for elem in discharge]
    EN = ['%.1f' % elem for elem in endarea]
    TO = ['%.1f' % elem for elem in topwidth]
    RS = ["%.4f" % reachslope]*len(ST)

    ratingdata = []
    for st, di, en, to, rs in zip(ST, DI, EN, TO, RS):
        ratingdata.append([st, di, en, to, rs])

    TWE_min = float(min(stage))
    TWE_max = float(max(stage))
    ratingtype = "XS"
    minstage = round(min(stage) + Dbf,2)
    reachno = str(rating)

    # ******************************************************************************************************
    # open "rating table output" file in text editor
    # [No. of Xsections = No. of Reaches]
    # ******************************************************************************************************

    CopyFeatures_management(os.path.join(optfolder, "transect.shp"), os.path.join(optfolder, "xpoint_reach" + str(rating) + ".shp"))
    CopyFeatures_management(os.path.join(optfolder, "line_trun.shp"), os.path.join(optfolder, "xline_reach" + str(rating) + ".shp"))

    SetParameterAsText(set_i, True); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps(plot_data)); set_i = set_i + 1
    SetParameterAsText(set_i, round(TWE_max,2)); set_i = set_i + 1
    SetParameterAsText(set_i, round(TWE_min,2)); set_i = set_i + 1
    SetParameterAsText(set_i, areami2_usda,); set_i = set_i + 1
    SetParameterAsText(set_i, reachslope); set_i = set_i + 1
    SetParameterAsText(set_i, round(Wbf,2)); set_i = set_i + 1
    SetParameterAsText(set_i, round(Dbf,2)); set_i = set_i + 1
    SetParameterAsText(set_i, minstage); set_i = set_i + 1
    SetParameterAsText(set_i, ratingtype); set_i = set_i + 1
    SetParameterAsText(set_i, reachno); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps(ratingdata)); set_i = set_i + 1



