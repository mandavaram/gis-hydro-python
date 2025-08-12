# -*- coding: utf-8 -*-
"""
Modified 05/2024

@author: Javier.Mardones
"""
from arcpy import (RasterToPolygon_conversion,
                   Intersect_analysis,
                   env,
                   Dissolve_management,
                   CalculateField_management,
                   GetParameterAsText,
                   Raster,
                   SearchCursor,
                   GetRasterProperties_management,
                   Describe,
                   Clip_analysis,
                   BuildRasterAttributeTable_management,
                   DeleteField_management,
                   AddField_management,
                   Buffer_analysis,
                   ListFields,
                   Shift_management,
                   Delete_management,
                   JoinField_management,
                   MosaicToNewRaster_management,
                   CreateFeatureclass_management,
                   Project_management,
                   SpatialReference,
                   SetParameterAsText)
from arcpy.sa import (Times,
                      Divide,
                      Log2,
                      Lookup,
                      FlowLength,
                      Int,
                      Con,
                      Combine,
                      Expand,
                      IsNull,
                      Minus,
                      SetNull,
                      GreaterThanEqual,
                      Extent,
                      Reclassify,
                      RemapValue,
                      ExtractByMask,
                      ZonalStatisticsAsTable)
from arcpy.da import SearchCursor as SearchCursorda
from arcpy.da import UpdateCursor as UpdateCursorda
from arcpy.management import GetCellValue
import os
import json

### Input

projectname = GetParameterAsText(0)
dem_layer = GetParameterAsText(1)
land_layer = GetParameterAsText(2)
soil_layer = GetParameterAsText(3)
hyd_cond = GetParameterAsText(4)
acc_thr = float(GetParameterAsText(5))
userstream = GetParameterAsText(6) == "true"

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])
directorygdb = os.path.join(directory, "data","gishydro.gdb")
demgdb = os.path.join(directory, "data","dem.gdb")
soilsgdb = os.path.join(directory, "data","soils.gdb")

noaaprecipgdb = os.path.join(directory, "data","noaaprecip.gdb")

prov = os.path.join(directorygdb, "md_phyisio_provinces")
province = os.path.join(directorygdb, "md_physio_boundaries")
mapstpm = os.path.join(directorygdb, "mapstpm")
prec_grid = os.path.join(noaaprecipgdb, "p2yr24ha")
limestonem = os.path.join(directorygdb, "limestonem")

optfolder = os.path.join(directory, "projects")
optfolder = os.path.join(optfolder, projectname)

sr_map = SpatialReference(4326)
sr_md = SpatialReference(26985)

# paths

dem = os.path.join(demgdb, dem_layer)
soils = os.path.join(soilsgdb, soil_layer)
landuse = os.path.join(optfolder, "landuse.tif")

rast = Raster(dem)
cellsize = int(rast.meanCellWidth)

# Environmental variables

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = os.path.join(demgdb, "neddem30")

# User streams (if selected) pre processing

if userstream:
    dem = os.path.join(optfolder, "dem_ext.tif")
    flowacc = os.path.join(optfolder, "flowacc_ext.tif")
    flowdir = os.path.join(optfolder, "flowdir_ext.tif")
    infstreams = os.path.join(optfolder, "infstr_ext.tif")

else:
    flowacc = os.path.join(directory, "data", "hydro{}.gdb".format(cellsize), "flowacc{}".format(cellsize))
    flowdir = os.path.join(directory, "data", "hydro{}.gdb".format(cellsize), "flowdir{}".format(cellsize))
    infstreams = os.path.join(directory, "data", "hydro{}.gdb".format(cellsize), "infstreams{}".format(cellsize))
                 

if not userstream:
    Clip_analysis(os.path.join(directorygdb, "nhd_streams_mr"), os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "nhd_strs.shp"))

# *******************************************************************************************************
# gage checking and selection
# *******************************************************************************************************

usgsgages = os.path.join(directorygdb, "usgsgagesm")
mdgages = os.path.join(directorygdb, "mdgagedstreams2024")
outletpoint = os.path.join(optfolder, "usr_pour_point.shp")

# added on 10-13-2017: Only using intersect tool to identify gauges which overlap watershed and mdgagugedstreams2016 file MODIFIED 5/2020
Intersect_analysis([mdgages,os.path.join(optfolder, "wshed.shp")], os.path.join(optfolder, "mask_ints.shp"),"ALL","#","INPUT")
JoinField_management(os.path.join(optfolder, "mask_ints.shp"),"GAGE_ID",usgsgages,"GAGE_ID","GAGE_ID")
Intersect_analysis([os.path.join(optfolder, "mask_ints.shp"),outletpoint], os.path.join(optfolder, "gauge_outlet.shp"),"ALL","#","INPUT")

gagelist = []
gagevalue = SearchCursor(os.path.join(optfolder, "gauge_outlet.shp"),"","","GAGE_ID","")
for g in gagevalue:
    gid = g.getValue("GAGE_ID")
    gagelist.append(gid)

# Creates layers for delineating subwatersheds in the future
CreateFeatureclass_management(optfolder, "addasstreams.shp", "POLYLINE", "", "ENABLED", "DISABLED", sr_md)
CreateFeatureclass_management(optfolder, "addasoutlets.shp", "POINT", "", "ENABLED", "DISABLED", sr_md)
CreateFeatureclass_management(optfolder, "addasreservoir.shp", "POINT", "", "ENABLED", "DISABLED", sr_md)

### DATA MANAGEMENT

basingrid = Raster(os.path.join(optfolder, "basingrid.tif"))
extentdist = cellsize*5
env.extent = Extent(basingrid.extent.XMin-extentdist, basingrid.extent.YMin-extentdist, basingrid.extent.XMax+extentdist, basingrid.extent.YMax+extentdist)

basinexpand = Expand(basingrid, 5, [1])

accum = Times(flowacc, basinexpand)
accum.save(os.path.join(optfolder, "flowacc.tif"))

shedtab = SearchCursor(basingrid, "", "", "Count", "")
for row in shedtab:
    basinarea = row.getValue("Count")

## project data for leaflet

areami2 = float((basinarea * pow(cellsize, 2)) / 2588881)  # conversion into sq miles
if acc_thr/640 > areami2:
    acc_thr = areami2*9/10

srcpixel = acc_thr/ pow(cellsize, 2) /0.000247
infstr_local = Con(accum >= int(srcpixel), accum)
infstr_local = Int(infstreams)
infstr_local.save(os.path.join(optfolder, "infstreams.tif"))

streams = GreaterThanEqual(os.path.join(optfolder, "infstreams.tif"), 1)
RasterToPolygon_conversion(streams, os.path.join(optfolder, "infstr_aux.shp"),"NO_SIMPLIFY","VALUE")
Dissolve_management(os.path.join(optfolder, "infstr_aux.shp"), os.path.join(optfolder, "infstr.shp"), "gridcode")
Project_management(os.path.join(optfolder, "infstr.shp"), os.path.join(optfolder, "infstr_proj.shp"), sr_map)

dem_clip = Times(dem, basinexpand)
dem_clip.save(os.path.join(optfolder, "dem_clip.tif"))

fdir = Times(flowdir, basinexpand)
fdir.save(os.path.join(optfolder, "flowdir.tif"))

soils_clip = Times(soils, basingrid)
soils_clip.save(os.path.join(optfolder, "soils.tif"))

## CURVE NUMBER RASTER

def get_lookup_file(lookup_type, condition):
    the_file = "{}lookup{}.txt".format(lookup_type, condition)
    return os.path.join(directory, "data", "lookup", the_file.lower())

lookup_dict = {
    "nlcd": {"nlcd2011", "nlcd2006", "nlcd2016", "nlcd2001", "nlcd2019"},
    "and": {"lu97m", "mdplu2002", "mdp2010"},
    "mdde": {"mdde2002"},
    "zoning": {"luult"},
    "mrlc": {"mrlc"},
    "usgs": {"lu70"}
}

for lookup_type, land_layers in lookup_dict.items():
    if land_layer in land_layers:
        lut_file = get_lookup_file(lookup_type, hyd_cond)
        break


lut_row = []
with open(lut_file, "r") as f:
    next(f)
    for line in f:
        lut_row.append(line.split("\t"))
lu_codes = [item[0] for item in lut_row]

mergerasts = Combine([landuse, os.path.join(optfolder, "soils.tif")])
mergerasts.save(os.path.join(optfolder,"cn_aux.tif"))

AddField_management(os.path.join(optfolder,"cn_aux.tif"), "cn", "LONG")
with UpdateCursorda(os.path.join(optfolder,"cn_aux.tif"),["landuse","soils","cn"]) as rows:
    for row in rows:
        if row[1] == -1:
            cnvalue = 100
        else:
            soilopts = lut_row[lu_codes.index(str(row[0]))][2:6]
            cnvalue = soilopts[int(row[1])-1]
        row[2] = str(cnvalue)
        rows.updateRow(row)
del rows

outRaster = Lookup(os.path.join(optfolder,"cn_aux.tif"),"cn")
outRaster.save(os.path.join(optfolder, "curvenumber.tif"))

try:
    Delete_management(os.path.join(optfolder, "cn_aux.tif"))
    Delete_management(os.path.join(optfolder, "infstr.shp"))
    Delete_management(os.path.join(optfolder, "wshed_aux.shp"))
    Delete_management(os.path.join(optfolder, "infstr_aux.shp"))
    Delete_management(os.path.join(optfolder, "mask_ints.shp"))
except:
    pass

### BASIN STATS

elevgrid = os.path.join(optfolder, "dem_clip.tif")
soils_path = os.path.join(optfolder, "soils.tif")
wats_lime = os.path.join(optfolder, "wats_lime.shp")
lime_int = os.path.join(optfolder, "lime_int.shp")
dirgrid = os.path.join(optfolder, "flowdir.tif")
landslope = os.path.join(optfolder, "landslope.tif")

### BASIN COMPOSITION

## convert all clipped rasters to polygons
lu_poly = RasterToPolygon_conversion(landuse, os.path.join(optfolder, "lu_poly.shp"), "NO_SIMPLIFY", "VALUE")
soil_poly = RasterToPolygon_conversion(soils_path, os.path.join(optfolder, "soil_poly.shp"), "NO_SIMPLIFY", "VALUE")

## intersect land use and soil to prepare two polygons: "lu_soil" and "lu_cn"
Intersect_analysis([lu_poly, soil_poly], os.path.join(optfolder, "lu_soil.shp"), "ALL", "#", "INPUT")

## dissolve above intersected polygon
Dissolve_management(os.path.join(optfolder, "lu_soil.shp"), os.path.join(optfolder, "lu_soil_diss.shp"), "GRIDCODE;GRIDCODE_1", "#", "MULTI_PART", "DISSOLVE_LINES")

## add filed to both of above dissolved polygons and compute area in acres
if not len(ListFields(os.path.join(optfolder, "lu_soil_diss.shp"), "area")) > 0:
    AddField_management(os.path.join(optfolder, "lu_soil_diss.shp"), "area", "FLOAT", 15, 4)
CalculateField_management(os.path.join(optfolder, "lu_soil_diss.shp"), "area", "!shape.area@acres!", "PYTHON")

# prepre a list of lu codes to feed into lu_description function in order to obtain matching descriptions list
lu_match = []
sc = SearchCursor(landuse, "", "", "VALUE", "")
for i in sc:
    v = i.getValue("VALUE")
    lu_match.append(v)

# create list of lists with zeroes
soil_acre_lists = [[0, 0, 0, 0] for i in range(len(lu_match))]

# preapre a list of soil acreage using lu_match list
lc_soil_diss = []
soil_lc_diss = []
lc_soil_aa = []
lu_soil_sc = SearchCursor(os.path.join(optfolder, "lu_soil_diss.shp"), "", "", "GRIDCODE;GRIDCODE_1;area", "")
for s in lu_soil_sc:
    lc = s.getValue("GRIDCODE")
    lc_soil_diss.append(lc)
    sc = s.getValue("GRIDCODE_1")
    soil_lc_diss.append(sc)
    aa = s.getValue("area")
    lc_soil_aa.append(round(aa, 2))

for idx, lu in enumerate(lu_match):
    for l, s, a in zip(lc_soil_diss, soil_lc_diss, lc_soil_aa):
        if l == lu:
            soil_acre_lists[idx][int(s) - 1] = a

# run hydro function to obtain land use description of categories present in watershed MODIFIED 2020
lu_code = []
lu_desc = []
with open(lut_file, "r") as f:
    next(f)
    for line in f:
        lu_strip = line.split("\t")[0] # lu code
        lu_code.append(lu_strip)
        dc_strip = line.split("\t")[1] # lu description
        lu_desc.append(dc_strip)
#convert lu_code to integer
lu_code = [int(i) for i in lu_code]
d_new = dict(zip(lu_code,lu_desc))
rev1 = {v:k for v,k in d_new.items()}
lu_desc = [rev1.get(item,item)  for item in lu_match]


# sum list of lists separately and cat at the end of lu description
total_area = [round(sum(i),2) for i in zip(*soil_acre_lists)]

# loop over land use, related total acreage, percent of land covered by this lu category, and A-B-C-D curve numbers
curve_num = []
for l in lu_match:
    with open(lut_file, "r") as f:
        next(f)
        for line in f:
            fields = line.split("\t")
            luc = fields[0]
            if int(l) == int(luc):
                # CN A, B, C, D are at indices 2, 3, 4, 5 respectively
                temp = [fields[i] for i in range(2, 6)]
                curve_num.append(temp)

# sum areas for each sub-list individually
acres = [round(sum(i),2) for i in soil_acre_lists]
total_all = sum(total_area)
area_percent = [round(float(ac_x / total_all) * 100, 2) for ac_x in acres]

sumcn = [sum(float(a)*float(b) for a,b in zip(*rows)) for rows in zip(soil_acre_lists, curve_num)]
avgCN = round(sum([float(b)*float(a)/100/float(c) for a,b,c in zip(area_percent,sumcn,acres)]),2)


### BASIN STATISTICS

env.snapRaster = os.path.join(demgdb, "neddem30")

# *******************************************************************************************************
# Warning messages
# *******************************************************************************************************
Impwarntext = """
        IMPERVIOUS AREA IN WATERSHED EXCEEDS 10%.
        Calculated discharges from USGS Regression
        Equations may not be appropriate.
                 """
provwarntext = """
        Watershed is within 5km of physiographic
        province boundary.  You should consider
        sensitivity of discharges to region location.
                 """
limewarntext = """
        Watershed is within 1km of underlying limestone
        geology.  You should consider sensitivity
        of discharges to percent limestone calculated.
                 """

# *******************************************************************************************************
# Get outlet coordinates and prepare masked grids for calculations
# *******************************************************************************************************
out_rast = Raster(elevgrid)
cellsize = out_rast.meanCellWidth
cellsq = cellsize * cellsize

# Get basingrid count [number of pixels]
shedtab = SearchCursor(basingrid, "", "", "Count", "")
for row in shedtab:
    basinarea = row.getValue("Count")

# *******************************************************************************************************
# Compute channel and land slope
# *******************************************************************************************************

# Multiply elevation grid and basin grid
elevation_clipped = elevgrid * basingrid

# Multiply direction grid and basin grid
basin_direction_grid = dirgrid * basingrid

# Compute downstream and upstream flow lengths
downstream_grid = FlowLength(basin_direction_grid, "DOWNSTREAM", "")
upstream_grid = FlowLength(basin_direction_grid, "UPSTREAM", "")

# Compute sum of downstream and upstream grids
sum_grid = downstream_grid + upstream_grid
sum_grid.save(os.path.join(optfolder, "maxlength.tif"))

# Get maximum value from upstream grid
max_length = GetRasterProperties_management(upstream_grid, "MAXIMUM")
max_length = max_length.getOutput(0)

# Compute tolerance
tolerance = 0.1 * cellsize

# Compute longest path
longest_path = Con(sum_grid > (float(max_length) - tolerance), 1 - sum_grid, 0)
longest_path = IsNull(longest_path)
longest_path = Con(longest_path == 0, 1, 0)

# Compute longest path upstream grid
longest_path_upstream_grid = longest_path * upstream_grid

# Compute thresholds for longest path
max_length_90_percent = 0.9 * float(max_length)
max_length_15_percent = 0.15 * float(max_length)

# Apply thresholds to longest path
longest_path_90_percent = Con(longest_path_upstream_grid < max_length_90_percent, longest_path_upstream_grid, 0)
longest_path_15_to_90_percent = Con(longest_path_90_percent > max_length_15_percent, longest_path_90_percent, 0)

# Set null values for longest path between 15% and 90%
longest_path_15_to_90_percent = SetNull(longest_path_15_to_90_percent, 1, "VALUE = 0")

# Multiply longest path between 15% and 90% with clipped elevation
longest_path_elevation = longest_path_15_to_90_percent * elevation_clipped

# Convert maximum length to miles
max_length_miles = float(max_length) * 0.000621371

# Get minimum and maximum elevation from longest path elevation
min_elevation = GetRasterProperties_management(longest_path_elevation, "MINIMUM").getOutput(0)
max_elevation = GetRasterProperties_management(longest_path_elevation, "MAXIMUM").getOutput(0)

# Compute slope in feet per mile
slope_feet_per_mile = (float(max_elevation) - float(min_elevation)) / (max_length_miles * 0.75)

# Convert slope to feet per foot
slope_feet_per_foot = slope_feet_per_mile / 5280.0


# *******************************************************************************************************
# Get LU count for Urban, Nil, Forest, and Storage -- it will be different for
# MOP, MD/DE, MRLC, and USGS
# *******************************************************************************************************
# *******************************************************************************************************
# Get Impervious count -- count will vary depending upon choice of input Landuse and Hyd condition
# *******************************************************************************************************

#*******************************************************************************************************
# Calculate LU cell count for delineated watershed -- LU count will vary based on landuse type
#*******************************************************************************************************

#*******************************************************************************************************
# Reclassify LU data to corresponding Imp values from NLCD table for Impervious area calculation
#*******************************************************************************************************

def get_lucount_file(lookup_type):
    the_file = "{}lucount.txt".format(lookup_type)
    return os.path.join(directory, "data", "lookup", the_file.lower())

def get_impcount_file(lookup_type):
    the_file = "{}impcount.txt".format(lookup_type)
    return os.path.join(directory, "data", "lookup", the_file.lower())

lookup_dict = {
    "nlcd": {"nlcd2011", "nlcd2006", "nlcd2016", "nlcd2001", "nlcd2019"},
    "and": {"lu97m", "mdplu2002", "mdp2010"},
    "mdde": {"mdde2002", "luult"},
    "mrlc": {"mrlc"},
    "usgs": {"lu70"}
}

for lookup_type, land_layers in lookup_dict.items():
    if land_layer in land_layers:
        lucount_file = get_lucount_file(lookup_type)
        impcount_file = get_impcount_file(lookup_type)
        break

lucount_rows = []
with open(lucount_file, "r") as f:
    next(f)
    for line in f:
        lucount_rows.append(line.split("\t"))

impcount_rows = []
with open(impcount_file, "r") as f:
    next(f)
    for line in f:
        impcount_rows.append(line.split("\t"))

# Initialize variables
UrbPct, FC, ST, IA = 0, 0, 0, 0

# Create a dictionary for faster lookup
lu_dict = {str(item[0]): item[1] for item in lucount_rows}
imp_dict = {str(item[0]): item[1] for item in impcount_rows}

# Iterate over rows
rows = SearchCursor(landuse, "", "", "Value;Count", "")
for row in rows:
    value = str(row.getValue("Value"))
    count = int(row.getValue("Count"))

    # Update variables based on category value
    catvalue = lu_dict.get(value)
    if catvalue == "1":
        UrbPct += count / basinarea * 100
    elif catvalue == "3":
        FC += count / basinarea * 100
    elif catvalue == "4":
        ST += count / basinarea * 100

    # Update IA based on impervious value
    impvalue = imp_dict.get(value)
    if impvalue is not None:
        IA += int(impvalue) * count / basinarea

del rows


# *******************************************************************************************************
# Get Limestone percent count
# *******************************************************************************************************

LIcnt = 0

#get geometry of single polygon in shapefile
with SearchCursorda(limestonem,['SHAPE@']) as cursor:
    for row in cursor: AOI_geom = row[0]

raster_extent = Describe(basingrid).extent
if raster_extent.overlaps(AOI_geom):

    limegrid = ExtractByMask(basingrid,limestonem)
    limegrid.save(os.path.join(optfolder, "limegrid.tif"))
    BuildRasterAttributeTable_management(os.path.join(optfolder, "limegrid.tif"), "Overwrite")
    with SearchCursorda(os.path.join(optfolder, "limegrid.tif"), "Count") as rows:
        for row in rows:
            LIcnt += row[0] or 0

LI = round(float((float(LIcnt) / basinarea) * 100),1)

areami2 = float((basinarea * cellsq) / 2588881)  # conversion into sq miles

# *******************************************************************************************************
# Get basin relief [it is difference of mean elevation and outlet elevation]
# *******************************************************************************************************
elev1 = GetRasterProperties_management(elevation_clipped, "MEAN")
mean_elev = float(elev1.getOutput(0))

for row in SearchCursorda(os.path.join(optfolder, "usr_pour_point.shp"), ["SHAPE@XY"]):
    x, y = row[0]
outlet_elev = GetCellValue(elevgrid, "{} {}".format(x, y))
outletelev = float(outlet_elev.getOutput(0))
basinrelief = float(mean_elev - outletelev)  # Assuming it is already converted into feets

# *******************************************************************************************************
# Get percent soil types from soil dataset
# *******************************************************************************************************

#*******************************************************************************************************
# Calculate soil percentages using selected soil data for soil percentages
#*******************************************************************************************************

def SoilPct(inSoil):
    soil_counts = {1: 0, 2: 0, 3: 0, 4: 0, -1: 0}
    temptab = SearchCursor(inSoil, "", "", "Value;Count", "")
    for row in temptab:
        thecount = float(row.getValue("Count"))
        value = row.getValue("Value")
        if value in soil_counts:
            soil_counts[value] += thecount

    return tuple(soil_counts.values())

pctR = SoilPct(soils_path)
pctSoil = list(map(float, pctR))
pctSoilR = [float((pct / basinarea) * 100) for pct in pctSoil]
pctAsoilR, pctBsoilR, pctCsoilR, pctDsoilR, pctWsoilR = pctSoilR


"""
The following code calculates the Time of Concentration. If multiple provinces
are involved, tc is weighted average of area of watershed in each province.
More correct would be to perform weighted average based on length of channel
in each province.  This modification will be performed at a later time.  (GEM - 12/01/99)
"""

# *******************************************************************************************************
# don"t add "theVTab" to TOC -- it will change list by drawing order to list
# by source which will prohibit addition of new layers
# *******************************************************************************************************
ZonalStatisticsAsTable(prov, "PROVINCE", basingrid, "theVTab", "DATA", "ALL")
DeleteField_management("theVTab", "ZONE_CODE;MIN;MAX;RANGE;MEAN;STD;SUM;VARIETY;MAJORITY;MINORITY;MEDIAN")
addFieldNameList = ["Q1.25", "Q1.50", "Q1.75", "Q2", "Q5", "Q10", "Q25", "Q50", "Q100", "Q200", "Q500"]
for each in addFieldNameList:
    if not len(ListFields("theVTab", each)) > 0:
        AddField_management("theVTab", each, "FLOAT", 10, 3)

# *******************************************************************************************************
sumarea = 0
theVTab = SearchCursor("theVTab", "", "", "Count", "")
for each in theVTab:
    count = each.getValue("Count")
    sumarea = sumarea + count
sumArea = sumarea
del each

# *******************************************************************************************************
# Compute Time of Concentration:
#                               1]  W.O. Thomas, Jr. Equation   [tc]
#                               2]  SCS Lag equation * 1.67     [lagtime]
# *******************************************************************************************************
sumtc = 0
theVTab = SearchCursor("theVTab", "", "", "Province;Count", "")
for row in theVTab:
    theProv = row.getValue("Province")
    theArea = row.getValue("Count")

    if row.getValue("Province") == "A":
        temptc = 0.133 * ((max_length_miles) ** (0.475)) * ((slope_feet_per_mile) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.194))
    elif row.getValue("Province") == "W":
        temptc = 0.133 * ((max_length_miles) ** (0.475)) * ((slope_feet_per_mile) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.366))
    elif row.getValue("Province") == "E":
        temptc = 0.133 * ((max_length_miles) ** (0.475)) * ((slope_feet_per_mile) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154)) * ((10) ** (0.366))
    else:
        temptc = 0.133 * ((max_length_miles) ** (0.475)) * ((slope_feet_per_mile) ** (-0.187)) * ((101 - FC) ** (-0.144)) * (
                (101 - IA) ** (0.861)) * ((ST + 1) ** (0.154))
    sumtc = sumtc + (temptc * theArea)
del row

tc = (sumtc / basinarea)

# *******************************************************************************************************
# Calculate landslope
# *******************************************************************************************************

def gishydroslope(eg, dg, cs, of, ext):
    dlgrid = Times(Con(((Log2(dg)) % 2) > 0, pow(2, 0.5), 1), cs)
    shifts = [(-cs, 0), (-cs, cs), (0, cs), (cs, cs), (cs, 0), (cs, -cs), (0, -cs), (-cs, -cs)]
    values = [1, 2, 4, 8, 16, 32, 64, 128]
    slope = 0

    for shift, value in zip(shifts, values):
        shift_temp = Divide(Minus(eg, Shift_management(eg, os.path.join(of, "dir_shift_{}.tif".format(value)), *shift)), dlgrid)
        slope = Con(dg, shift_temp, slope, "VALUE={}".format(value))

    # Slope value was going below 0 which isn't correct (raw dem used instead of filled). Condition is added to at least have 0.0001
    slopegrid = Con(slope > 0.0001, slope, 0.0001)
    slopegrid = Times(slopegrid, basingrid)
    slopegrid.save(os.path.join(of, "landslope" + ext + ".tif"))
    return slopegrid


# *******************************************************************************************************
# Calculate lagtime
# *******************************************************************************************************
landslopegrid = gishydroslope(elevgrid,dirgrid,cellsize,optfolder, "")

landsloperesult = GetRasterProperties_management(landslopegrid, "MEAN")
landslopevalue = float(landsloperesult.getOutput(0))
landslope = float(landslopevalue) / 3.28084  # modified: 09/06/2018 (divided by 3.28084)

lagtime = (float(100 * ((max_length_miles * 5280) ** (0.8)) * (((1000 / avgCN) - 9) ** (0.7))) / float(1900 * ((abs(landslope) * 100) ** (0.5)))) / 60

# *******************************************************************************************************
# Calculate Mean Annual Precipitation
# *******************************************************************************************************

env.snapRaster = os.path.join(demgdb, "neddem30")
env.cellSize = elevgrid

maprecbasin = Times(mapstpm, basingrid)  # Make sure basingrid has value 1 otherwise all precip will be 0
theprec = Times(prec_grid, basingrid)
precavg = GetRasterProperties_management(maprecbasin, "MEAN")
precavg = float(precavg.getOutput(0))
avgprec = GetRasterProperties_management(theprec, "MEAN")
avgprec = float(avgprec.getOutput(0))
maprec = float(precavg / (1000 * 2.54))
p2yr = float(avgprec / 1000)


# *******************************************************************************************************
# Print out Impervious area warning message
# *** warning message is included in for loop despite the fact that technically it could be printed
#     twice. Since both "Appalachian Plateau" and "Eastern Coastal Plain" are far apart so it is
#     impossible to have that big watershed while doing analysis with GISHydroNXT
# *******************************************************************************************************
html_warning = ""
theVTab = SearchCursor("theVTab", "", "", "Province", "")
for row in theVTab:
    if (row.getValue("Province") == "A") or (row.getValue("Province") == "E"):
        if float(IA) >= 10:
            html_warning = html_warning + Impwarntext + "\n"

# *******************************************************************************************************
# Close to boundary condition for provinces -- Near tool isn"t available with
# basic level license therefore a more crude method was emplyed here. It can
# be improved in future by using "Geometry()" tool to get distance
# *******************************************************************************************************
Buffer_analysis(os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "wats_prov.shp"), "5000", "#", "#", "ALL", "FID")
Intersect_analysis([province, os.path.join(optfolder, "wats_prov.shp")], os.path.join(optfolder, "prov_int.shp"))
prov_cursor = SearchCursor(os.path.join(optfolder, "prov_int.shp"), "", "", "FID", "")
prov = prov_cursor.next()
if prov != None:
    html_warning = html_warning + provwarntext + "\n"

# *******************************************************************************************************
# Close to boundary condition for limestone -- Near tool isn"t available with
# basic level license therefore a more crude method was emplyed here. It can
# be improved in future by using "Geometry()" tool to get distance
# *******************************************************************************************************
Buffer_analysis(os.path.join(optfolder, "wshed.shp"), os.path.join(optfolder, "wats_lime.shp"), "1000", "#", "#", "ALL", "FID")
Intersect_analysis([limestonem, wats_lime], lime_int, "ALL", "#", "INPUT")
lime_cursor = SearchCursor(lime_int, "", "", "FID", "")
lime = lime_cursor.next()
if lime != None:
    html_warning = html_warning + limewarntext + "\n"

# *******************************************************************************************************
# Hydraulic coefficients for channel geometry -- default to outlet"s physiographic province
# *******************************************************************************************************

coefficients = {
    "A": {"a": 13.87, "b": 0.44, "c": 0.95, "d": 0.31, "e": 13.17, "f": 0.75},
    "B": {"a": 13.87, "b": 0.44, "c": 0.95, "d": 0.31, "e": 13.17, "f": 0.75},
    "P": {"a": 14.78, "b": 0.39, "c": 1.18, "d": 0.34, "e": 17.42, "f": 0.73},
    "W": {"a": 10.3, "b": 0.38, "c": 1.01, "d": 0.32, "e": 10.34, "f": 0.70},
    "E": {"a": 10.3, "b": 0.38, "c": 1.01, "d": 0.32, "e": 10.34, "f": 0.70}
}

province_names = {
    "A": "Appalachian Plateaus and Allegheny Ridges",
    "B": "Blue Ridge and Piedmont",
    "P": "Blue Ridge and Piedmont",
    "W": "Western Coastal Plain",
    "E": "Eastern Coastal Plain"
}

provstring = []
theVTab = SearchCursor("theVTab", "", "", "Province;Count", "")

for row in theVTab:
    province = row.getValue("Province")
    area_field = float(row.getValue("Count"))
    prov_percent = (area_field / sumArea) * 100

    if province in coefficients:
        coef = coefficients[province]
        Coef_W, Exp_W, Coef_D, Exp_D, Coef_A, Exp_A = coef.values()

    if province in province_names:
        provstring.append([province_names[province], "{0:.2f}".format(prov_percent)])
    else:
        provstring.append(["No Province Selected", "0"])

coef_list = [Coef_W, Coef_D, Coef_A]
exp_list = [Exp_W, Exp_D, Exp_A]

pctsoilR = [round(pctAsoilR,2), round(pctBsoilR,2), round(pctCsoilR,2), round(pctDsoilR,2)]

# *******************************************************************************************************
# format precision before text file string settings
# *******************************************************************************************************
max_length_miles = "{0:.2f}".format(max_length_miles)
slope_feet_per_mile = "{0:.2f}".format(slope_feet_per_mile)
slope_feet_per_foot = "{0:.8f}".format(slope_feet_per_foot)
landslope = "{0:.8f}".format(landslope)
UrbPct = "{0:.2f}".format(UrbPct)
FC = "{0:.2f}".format(FC)
ST = "{0:.2f}".format(ST)
IA = "{0:.2f}".format(IA)
LI = "{0:.2f}".format(LI)
areami2 = "{0:.2f}".format(areami2)
outletelev = "{0:.2f}".format(outletelev)
basinrelief = "{0:.2f}".format(basinrelief)
avgCN = "{0:.1f}".format(avgCN)
pctAsoilR = "{0:.2f}".format(pctAsoilR)
pctBsoilR = "{0:.2f}".format(pctBsoilR)
pctCsoilR = "{0:.2f}".format(pctCsoilR)
pctDsoilR = "{0:.2f}".format(pctDsoilR)
tc = "{0:.2f}".format(tc)
lagtime = "{0:.2f}".format(lagtime)
maprec = "{0:.2f}".format(maprec)
p2yr = "{0:.2f}".format(p2yr)


### FRRE DATA

landslope_frre = float(landslope)
if int(cellsize) == 10:
    elevgrid = Times(basingrid, os.path.join(directory, "data","dem.gdb", "neddem30"))
    dirgrid = Times(basingrid, os.path.join(directory, "data", "hydro30.gdb", "flowdir30"))
    
    landslopegrid = gishydroslope(elevgrid,dirgrid,cellsize, optfolder ,"_frre")
    landsloperesult = GetRasterProperties_management(landslopegrid, "MEAN")
    landslopevalue = float(landsloperesult.getOutput(0))
    landslope_frre = float(landslopevalue) / 3.28084


landusegdb = os.path.join(directory, "data","landuse.gdb")
landuse_mdp = os.path.join(landusegdb, "mdp2010")
landuse_nlcd = os.path.join(landusegdb, "nlcd2019")
        
env.extent = basingrid
env.snapRaster = os.path.join(demgdb, "neddem30")

landuse_mdp = Times(landuse_mdp, basingrid)
landuse_mdp.save(os.path.join(optfolder, "landuse_mdp.tif"))
landuse_mdp = os.path.join(optfolder, "landuse_mdp.tif")
BuildRasterAttributeTable_management(landuse_mdp, "Overwrite")

landuse_nlcd = Times(landuse_nlcd, basingrid)
landuse_nlcd.save(os.path.join(optfolder, "landuse_nlcd.tif"))
landuse_nlcd = os.path.join(optfolder, "landuse_nlcd.tif")
BuildRasterAttributeTable_management(landuse_nlcd, "Overwrite")

mdp_impcount_file = os.path.join(directory, "data","lookup","andimpcount.txt")
nlcd_impcount_file = os.path.join(directory, "data","lookup","nlcdimpcount.txt")
mdp_raster_loc = os.path.join(optfolder, "mdp_imp.tif")
nlcd_raster_loc = os.path.join(optfolder, "nlcd_imp.tif")

# Load the text file and create a dictionary mapping old values to new values
def create_impraster(impcount_file, landuseraster, raster_loc):
    value_map = []
    with open(impcount_file, 'r') as f:
        next(f)
        for line in f:
            old_value, new_value = line.strip().split('\t')
            value_map.append([int(old_value), int(new_value)])
    new_raster = Reclassify(Raster(landuseraster), "Value", RemapValue(value_map))
    new_raster.save(raster_loc)

create_impraster(nlcd_impcount_file, landuse_nlcd, nlcd_raster_loc)

try:
    create_impraster(mdp_impcount_file, landuse_mdp, mdp_raster_loc)    
    MosaicToNewRaster_management([Raster(mdp_raster_loc), Raster(nlcd_raster_loc)], optfolder, "imp_frre.tif", sr_md, "8_BIT_UNSIGNED", "30", 1, "FIRST")
except:
    MosaicToNewRaster_management([Raster(nlcd_raster_loc)], optfolder, "imp_frre.tif", sr_md, "8_BIT_UNSIGNED", "30", 1, "FIRST")

IA_frre = GetRasterProperties_management(os.path.join(optfolder, "imp_frre.tif"), "MEAN")
IA_frre = float(IA_frre.getOutput(0))
  
try:
    Delete_management(landuse_mdp)
    Delete_management(landuse_nlcd)
except:
    pass
  
#*******************************************************************************************************
# Calculate soil percentages using SSURGO for use in regression equations
#*******************************************************************************************************
# soil types are initialized with zero values before function arguments to avoid following error:
# UnboundLocalError: local variable referenced before assignment
#

ssurgo = Times(os.path.join(soilsgdb, "ssurgo_2021.tif"), basingrid)

pctAsoil = 0
temptab = SearchCursor(ssurgo,"","","Value;Count","")
for row in temptab:
    count = int(row.getValue("Count"))
    if row.getValue("Value") == 1:
        pctAsoil = pctAsoil + count

pctAsoil_frre = float(pctSoil[0]/basinarea*100)

frre_param = [round(landslope_frre,5), round(IA_frre,2), round(pctAsoil_frre,2)]

#############
## OUTPUTS ##
#############

# Data
SetParameterAsText(7, json.dumps(gagelist))
SetParameterAsText(8, os.path.join(optfolder, "infstr_proj.shp"))

# Basin comp
SetParameterAsText(9, json.dumps(lu_desc))
SetParameterAsText(10, json.dumps(soil_acre_lists))
SetParameterAsText(11, json.dumps(total_area))
SetParameterAsText(12, json.dumps(acres))
SetParameterAsText(13, json.dumps(area_percent))
SetParameterAsText(14, json.dumps(curve_num))

# Basins statistics
SetParameterAsText(15, html_warning)
SetParameterAsText(16, int(x))
SetParameterAsText(17, int(y))
SetParameterAsText(18, json.dumps(provstring))
SetParameterAsText(19, areami2)
SetParameterAsText(20, slope_feet_per_mile)
SetParameterAsText(21, slope_feet_per_foot)
SetParameterAsText(22, landslope)
SetParameterAsText(23, UrbPct)
SetParameterAsText(24, IA)
SetParameterAsText(25, tc)
SetParameterAsText(26, lagtime)
SetParameterAsText(27, max_length_miles)
SetParameterAsText(28, outletelev)
SetParameterAsText(29, basinrelief)
SetParameterAsText(30, avgCN)
SetParameterAsText(31, FC)
SetParameterAsText(32, ST)
SetParameterAsText(33, LI)
SetParameterAsText(34, json.dumps(pctsoilR))
SetParameterAsText(35, json.dumps(frre_param))
SetParameterAsText(36, p2yr)
SetParameterAsText(37, maprec)
SetParameterAsText(38, json.dumps(coef_list))
SetParameterAsText(39, json.dumps(exp_list))

try:

    Delete_management(os.path.join(optfolder, "lu_poly.shp"))
    Delete_management(os.path.join(optfolder, "lu_soil.shp"))
    Delete_management(os.path.join(optfolder, "lu_soil_diss.shp"))
    Delete_management(os.path.join(optfolder, "soil_poly.shp"))

    Delete_management(os.path.join(optfolder, "dir_shift_1.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_2.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_4.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_8.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_16.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_32.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_64.tif"))
    Delete_management(os.path.join(optfolder, "dir_shift_128.tif"))
    Delete_management(os.path.join(optfolder, "lime_int.shp"))
    Delete_management(os.path.join(optfolder, "prov_int.shp"))
    Delete_management(os.path.join(optfolder, "wats_lime.shp"))
    Delete_management(os.path.join(optfolder, "wats_prov.shp"))
    Delete_management(os.path.join(optfolder, "gauge_outlet.shp"))
except:
    pass

