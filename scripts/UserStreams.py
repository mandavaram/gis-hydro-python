# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 13:21:14 2022

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   Extent,
                   Array,
                   Point,
                   Polyline,
                   PointGeometry,
                   Polygon,
                   Clip_analysis,
                   PolylineToRaster_conversion,
                   Describe,
                   RasterToPolygon_conversion,
                   Merge_management,
                   CreateFeatureclass_management,
                   CopyFeatures_management,
                   SpatialReference,
                   Project_management,
                   Clip_management,
                   Dissolve_management,
                   SetParameterAsText)
from arcpy.sa import (IsNull,
                      Minus,
                      Times,
                      Plus,
                      Fill,
                      FlowDirection,
                      FlowAccumulation,
                      Con)
from arcpy.da import InsertCursor as InsertCursorda
import os
import json
import time

# Constants
DATA_DIR = "data"
SR_MAP_WKID = 4326
SR_MD_WKID = 26985

# Parameters
proj_name = GetParameterAsText(0).replace(" ", "_")
dem_layer = GetParameterAsText(1)
extent = GetParameterAsText(2)
userflag = GetParameterAsText(3) == "true"
usrstreams = GetParameterAsText(4)
nhdopt = GetParameterAsText(5)
acc_thr = float(GetParameterAsText(6))

# Paths
directory = os.path.join(os.environ['GISHydro_DIR'])
directorygdb = os.path.join(directory, DATA_DIR, "gishydro.gdb")
demgdb = os.path.join(directory, DATA_DIR, "dem.gdb")
dem = os.path.join(demgdb, dem_layer)

# Spatial references
sr_map = SpatialReference(SR_MAP_WKID)
sr_md = SpatialReference(SR_MD_WKID)

# Create temporary folder if it doesn't exist
temp_folder = os.path.join(directory, "projects")
os.makedirs(temp_folder, exist_ok=True)

# Create project folder
folder_name = f"{time.strftime('%Y%m%d_%H%M%S')}_{proj_name}"
project_folder = os.path.join(temp_folder, folder_name)
os.makedirs(project_folder, exist_ok=True)

# Set environment variables
env.overwriteOutput = True
env.scratchWorkspace = project_folder
env.workspace = project_folder
env.snapRaster = dem

# Function to transform coordinates
def coord_transf(x1, y1, sr1, sr2):
    pt1 = Point(x1, y1)
    ptgeo1 = PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    return float(pt2.X), float(pt2.Y)

# Save AOI extent
extent = json.loads(extent)
coordinates = extent['geometry']['coordinates'][0]
east, north = coord_transf(float(coordinates[0][0]), float(coordinates[1][1]), sr_map, sr_md)
west, south = coord_transf(float(coordinates[2][0]), float(coordinates[3][1]), sr_map, sr_md)

extent_clip = Extent(east, south, west, north)
env.extent = extent_clip

# Create mask layer
array = Array([Point(extent_clip.XMin, extent_clip.YMin),
               Point(extent_clip.XMin, extent_clip.YMax),
               Point(extent_clip.XMax, extent_clip.YMax),
               Point(extent_clip.XMax, extent_clip.YMin),
               Point(extent_clip.XMin, extent_clip.YMin)])
mask_layer = Polygon(array, sr_md)

# Clip DEM
dem_clip = os.path.join(project_folder, "dem_ext.tif")
mask = os.path.join(project_folder, "mask.shp")
CopyFeatures_management(mask_layer, mask)
Clip_management(dem, "", dem_clip, mask)

# Create user streams feature class
CreateFeatureclass_management(project_folder, "userstreams.shp", "POLYLINE", "", "ENABLED", "DISABLED", sr_md)

# Process user streams
if userflag:
    usrstreams = json.loads(usrstreams)
    for i, linefeatures in enumerate(usrstreams['features']):
        coordinates = linefeatures['geometry']['coordinates']
        arcpoints = [Point(*coord_transf(float(coord[0]), float(coord[1]), sr_map, sr_md)) for coord in coordinates]
        polyline = Polyline(Array(arcpoints))
        insertcursor = InsertCursorda(os.path.join(project_folder, "userstreams.shp"), ("OID@", "SHAPE@"))
        insertcursor.insertRow([i, polyline])

# Process NHD streams
if nhdopt in ["1", "2"]:
    nhd = os.path.join(directorygdb, "nhd_streams_hr" if nhdopt == "2" else "nhd_streams_mr")
    Clip_analysis(nhd, mask, os.path.join(project_folder, "nhd_strs.shp"))
    Merge_management([os.path.join(project_folder, "userstreams.shp"), os.path.join(project_folder, "nhd_strs.shp")], os.path.join(project_folder, "burn.shp"))
    try:
        PolylineToRaster_conversion(os.path.join(project_folder, "burn.shp"), "FID", os.path.join(project_folder, "nhd_rast.tif"), "", "", dem)
    except:
        PolylineToRaster_conversion(os.path.join(project_folder, "burn.shp"), "OBJECTID", os.path.join(project_folder, "nhd_rast.tif"), "", "", dem)
    burn = 100
    calc1 = IsNull(os.path.join(project_folder, "nhd_rast.tif"))
    calc2 = Minus(calc1, 1)
    calc1 = Times(calc2, -1 * burn)
    calc2 = Minus(dem_clip, calc1)
    burned = Plus(calc2, burn)
    burned.save(os.path.join(project_folder, "burned.tif"))
    fill = Fill(os.path.join(project_folder, "burned.tif"))
    fill.save(os.path.join(project_folder, "fill.tif"))
else:
    if userflag:
        try:
            PolylineToRaster_conversion(os.path.join(project_folder, "userstreams.shp"), "FID", os.path.join(project_folder, "user_rast.tif"), "", "", dem)
        except:
            PolylineToRaster_conversion(os.path.join(project_folder, "userstreams.shp"), "OBJECTID", os.path.join(project_folder, "user_rast.tif"), "", "", dem)
        burn = 100
        calc1 = IsNull(os.path.join(project_folder, "user_rast.tif"))
        calc2 = Minus(calc1, 1)
        calc1 = Times(calc2, -1 * burn)
        calc2 = Minus(dem_clip, calc1)
        burned = Plus(calc2, burn)
        burned.save(os.path.join(project_folder, "burned.tif"))
        fill = Fill(os.path.join(project_folder, "burned.tif"))
        fill.save(os.path.join(project_folder, "fill.tif"))
    else:
        fill = Fill(dem_clip)
        fill.save(os.path.join(project_folder, "fill.tif"))

# Calculate flow direction and accumulation
flowdir = FlowDirection(fill)
flowdir.save(os.path.join(project_folder, "flowdir_ext.tif"))
flowacc = FlowAccumulation(flowdir)
flowacc.save(os.path.join(project_folder, "flowacc_ext.tif"))

# Calculate source pixel
faccdesc = Describe(flowacc)
cellSize = faccdesc.meanCellHeight
srcpixel = acc_thr / (cellSize ** 2) / 0.000247
infstr = Con(flowacc >= int(srcpixel), flowacc)
streams = infstr >= 1
streams.save(os.path.join(project_folder, "infstr_ext.tif"))

# Convert raster to polygon
RasterToPolygon_conversion(streams, os.path.join(project_folder, "infstr_aux.shp"), "NO_SIMPLIFY", "VALUE")
Dissolve_management(os.path.join(project_folder, "infstr_aux.shp"), os.path.join(project_folder, "infstr.shp"), "gridcode")
Project_management(os.path.join(project_folder, "infstr.shp"), os.path.join(project_folder, "infstr_proj.shp"), sr_map)

# Set output parameters
SetParameterAsText(7, folder_name)
SetParameterAsText(8, os.path.join(project_folder, "infstr_proj.shp"))