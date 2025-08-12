# -*- coding: utf-8 -*-
"""
Modified 02/2021

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   Raster,
                   SpatialReference,
                   PointToRaster_conversion,
                   Merge_management,
                   FeatureToRaster_conversion,
                   RasterToPoint_conversion,
                   RasterToPolygon_conversion,
                   Dissolve_management,
                   FieldMappings,
                   SpatialJoin_analysis,
                   SearchCursor,
                   Describe,
                   CreateFeatureclass_management,
                   ListFields,
                   AddField_management,
                   UpdateCursor,
                   DeleteField_management,
                   PolygonToRaster_conversion,
                   Project_management,
                   Delete_management,
                   GetCount_management,
                   SetParameterAsText)
from arcpy.sa import (Plus,
                      Con,
                      SetNull,
                      IsNull,
                      Times,
                      StreamLink,
                      ZonalStatistics,
                      Watershed,
                      StreamToFeature)
from arcpy.da import SearchCursor as SearchCursorda
from arcpy.da import InsertCursor as InsertCursorda
import os

projectname = GetParameterAsText(0)

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

optfolder = os.path.join(directory, "projects", projectname)

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = os.path.join(optfolder, "dem_clip.tif")
env.extent = os.path.join(optfolder, "dem_clip.tif")

for row in SearchCursorda(os.path.join(optfolder, "usr_pour_point.shp"), ["SHAPE@XY"]):
    x, y = row[0]

rast = Raster(os.path.join(optfolder, "dem_clip.tif"))
cellsize = rast.meanCellWidth

try:
    Delete_management(os.path.join(optfolder, "outlets_user.tif"))
except:
    pass

count = int(GetCount_management(os.path.join(optfolder, "addasoutlets.shp")).getOutput(0))
if count > 0:
    PointToRaster_conversion(os.path.join(optfolder, "addasoutlets.shp"), "FID", os.path.join(optfolder, "Outlets_temp.tif"), "MOST_FREQUENT", "NONE", cellsize)
    outlets_adj = Plus(os.path.join(optfolder, "Outlets_temp.tif"), 1)
    outlets_adj.save(os.path.join(optfolder, "AddOutlets.tif"))
    outlets_custom = Con(IsNull(os.path.join(optfolder, "AddOutlets.tif")), 0, os.path.join(optfolder, "AddOutlets.tif"))
    outlets_custom = SetNull(outlets_custom, outlets_custom, "VALUE = 0")
    outlets_custom.save(os.path.join(optfolder, "outlets_user.tif"))

env.extent = os.path.join(optfolder, "dem_clip.tif")
Merge_management(os.path.join(optfolder, "AddasStreams.shp"), os.path.join(optfolder, "StrmMerge.shp"))  # no apparent benefit of using merge
FeatureToRaster_conversion(os.path.join(optfolder, "StrmMerge.shp"), "Id", os.path.join(optfolder, "ModStr.tif"), "#")
Streams = Con(Raster(os.path.join(optfolder, "ModStr.tif")) == 0, 1)
Streams.save(os.path.join(optfolder, "ModStreams.tif"))  # add it to TOC

env.extent = os.path.join(optfolder, "dem_clip.tif")
flowacc = os.path.join(optfolder, "flowacc.tif")
flwdir = os.path.join(optfolder, "flowdir.tif")
strlnk = StreamLink(os.path.join(optfolder, "ModStreams.tif"), flwdir)
zonemax = ZonalStatistics(strlnk, "Value", flowacc, "MAXIMUM", "NODATA")
faccmax = Con(flowacc == zonemax, zonemax, IsNull(zonemax))
outlets = Con(faccmax > 0, faccmax)
outlets.save(os.path.join(optfolder, "outlets.tif"))

# if outlets are added by user then add mod stream and user specified outlets
env.extent = os.path.join(optfolder, "dem_clip.tif")
outlets_user   = os.path.join(optfolder, "outlets_user.tif")
if os.path.exists(outlets_user):
    aux_out1 = RasterToPoint_conversion(outlets)
    aux_out2 = RasterToPoint_conversion(outlets_user)
    aux_out3 = Merge_management([aux_out1,aux_out2])
    outlets = PointToRaster_conversion(aux_out3, "FID", os.path.join(optfolder, "outlets.tif"), "#", "#", cellsize)
    del aux_out1, aux_out2, aux_out3

subshed = os.path.join(optfolder, "subshed.shp")

env.extent = os.path.join(optfolder, "dem_clip.tif")
subwshed = Watershed(flwdir, os.path.join(optfolder, "outlets.tif"), "VALUE")
subwshed = Times(subwshed, os.path.join(optfolder, "basingrid.tif"))
RasterToPolygon_conversion(subwshed, os.path.join(optfolder, "tmpsubwshd.shp"),"NO_SIMPLIFY","VALUE")

env.extent = os.path.join(optfolder, "dem_clip.tif")
tmpsub = os.path.join(optfolder, "tmpsubwshd.shp")
Dissolve_management(tmpsub,subshed,"GRIDCODE","#","MULTI_PART","DISSOLVE_LINES")

subrivers = os.path.join(optfolder, "subrivers.shp")
# new stream link raster to handle added outlets
if os.path.exists(outlets_user):
    StrmFDRgrid = Times(os.path.join(optfolder, "ModStreams.tif"),flwdir)
    newLnkGrid = Watershed(StrmFDRgrid, outlets, "VALUE")
    newLnkGrid = Times(os.path.join(optfolder, "ModStreams.tif"),newLnkGrid)
    StreamToFeature(newLnkGrid, flwdir, subrivers, "NO_SIMPLIFY")
    del newLnkGrid
else:
    StreamToFeature(os.path.join(optfolder, "ModStreams.tif"), flwdir, subrivers, "NO_SIMPLIFY")
del outlets,flwdir

# *******************************************************************************************************
# correct sub-basin indexing -- order of sub-basins has to be the same as
# appearing in "subriver.shp"
# *******************************************************************************************************

spatial_join1 = os.path.join(optfolder, "subriver_subshed.shp")

# add table using spatial join (target = subriver, join = subshed)
fieldmappings = FieldMappings()
fieldmappings.addTable(subrivers)
fieldmappings.addTable(subshed)
SpatialJoin_analysis(subrivers, subshed, spatial_join1, "#", "#", fieldmappings, "HAVE_THEIR_CENTER_IN", "#", "#")

# get subriver FID and centroid (or midpoint) XY values
sub_cur = SearchCursor(spatial_join1, "", "", "Shape;ARCID;GRIDCODE", "")
x = []
y = []
for row in sub_cur:
    XMidPoint = row.shape.positionAlongLine(0.50, True).firstPoint.X
    x.append(XMidPoint)
    YMidPoint = row.shape.positionAlongLine(0.50, True).firstPoint.Y
    y.append(YMidPoint)

del sub_cur

xy = list(zip(x, y))
# create a point shapefile, add fields "ARCID" & "GRIDCODE", and add midpoint XY values to field "SHAPE@"
spatial_reference = Describe(spatial_join1).spatialReference
CreateFeatureclass_management(optfolder, "subriver_xy.shp", "POINT", "", "ENABLED", "DISABLED", spatial_reference)
subriver_xy = os.path.join(optfolder, "subriver_xy.shp")

# added on 09-30-2014 to debug field type error -- changed field type to "Double" with precision "10"
if not len(ListFields(subriver_xy, "ARCID")) > 0:
    AddField_management(subriver_xy, "ARCID", "DOUBLE", 10, "")
    AddField_management(subriver_xy, "GRIDCODE", "DOUBLE", 10, "")
sr_xy = InsertCursorda(subriver_xy, ("SHAPE@XY"))
for i in xy:
    sr_xy.insertRow([i])

# update "ARCID" and "GRIDCODE" row values using "subriver_subshed" shapefile
sub_cur = SearchCursor(spatial_join1, "", "", "ARCID;GRIDCODE", "")
ic = UpdateCursor(subriver_xy, "", "", "ARCID;GRIDCODE", "")

row1 = sub_cur.next()
row2 = ic.next()

while row1:
    row2.setValue("ARCID", row1.getValue("ARCID"))
    ic.updateRow(row2)
    row2.setValue("GRIDCODE", row1.getValue("GRIDCODE"))
    ic.updateRow(row2)
    row1 = sub_cur.next()
    row2 = ic.next()
del sub_cur, ic

# perform spatial join (target = subshed, join = subriver_xy) and lists extraction
target2 = os.path.join(optfolder, "subshed.shp")
joint2 = os.path.join(optfolder, "subriver_xy.shp")
spatial_join2 = os.path.join(optfolder, "subshed_subriver.shp")

fieldmappings = FieldMappings()
fieldmappings.addTable(target2)
fieldmappings.addTable(joint2)
SpatialJoin_analysis(target2, joint2, spatial_join2, "#", "#", fieldmappings, "INTERSECT", "#", "#")

# create a polygon which will contain subriver FID, ARCID, and GRIDCODE (obtained from spatially joined shapefile)
spatial_reference = Describe(os.path.join(optfolder, "subshed_subriver.shp")).spatialReference
CreateFeatureclass_management(optfolder, "subshed.shp", "POLYGON", "", "DISABLED", "DISABLED", spatial_reference)

sc = SearchCursorda(spatial_join2, ("ARCID", "GRIDCODE", "SHAPE@"))
storage = []
for row in sc:  # changed sF to sc [could be a mistake to have "sF"]
    storage.append(row)

sortedlist = sorted(storage, key=lambda x: x[0])  # list sort is based on "ARCID"

poly = os.path.join(optfolder, "subshed.shp")

# added on 09-30-2014 to debug field type error -- changed field type to "Double" with precision "10"
if not len(ListFields(poly, "ARCID")) > 0:
    AddField_management(poly, "ARCID", "DOUBLE", 10, "")
    AddField_management(poly, "GRIDCODE", "DOUBLE", 10, "")

DeleteField_management(poly, "Id")
fcout = InsertCursorda(poly, ("ARCID", "GRIDCODE", "SHAPE@"))

for i in sortedlist:
    fcout.insertRow(i)
del fcout, sc

# convert newly created and indexed sub-watershed shapefile to raster
PolygonToRaster_conversion(subshed, "ARCID", os.path.join(optfolder, "subsheds.tif"), "", "ARCID", os.path.join(optfolder, "dem_clip.tif"))

sr_map = SpatialReference(4326)
Project_management(os.path.join(optfolder, "subshed.shp"), os.path.join(optfolder, "subshed_proj.shp"), sr_map)

try:
    Delete_management(os.path.join(optfolder, "Outlets_temp.tif"))
    Delete_management(os.path.join(optfolder, "ModStr.tif"))
    Delete_management(os.path.join(optfolder, "ModStreams.tif"))
    Delete_management(os.path.join(optfolder, "outlets_user.tif"))
    Delete_management(os.path.join(optfolder, "StrmMerge.shp"))
    Delete_management(os.path.join(optfolder, "subriver_xy.shp"))
    Delete_management(os.path.join(optfolder, "subshed_subriver.shp"))
    Delete_management(os.path.join(optfolder, "tmpsubwshd.shp"))
except:
    pass

SetParameterAsText(1, os.path.join(optfolder, "subshed_proj.shp"))