# -*- coding: utf-8 -*-
"""
Modified 05/2024

@author: Javier.Mardones
"""

import os
import time
from arcpy import (env,
                   RasterToPolygon_conversion,
                   SpatialReference,
                   BuildRasterAttributeTable_management,
                   management,
                   RasterToPoint_conversion,
                   Raster,
                   Dissolve_management,
                   JoinField_management,
                   Clip_analysis,
                   GetParameterAsText,
                   Point,
                   Project_management,
                   Delete_management,
                   PointGeometry,
                   Clip_management,
                   SearchCursor,
                   SetParameterAsText)
from arcpy.sa import (IsNull,
                      Times,
                      Watershed,
                      SnapPourPoint,
                      Extent,
                      Con)
from arcpy.da import SearchCursor as SearchCursorda

# Constants
DATA_DIR = "data"
HYDRO_GDB_TEMPLATE = "hydro{}.gdb"
SR_MAP_CODE = 4326
SR_MD_CODE = 26985

# Parameters
proj_name = GetParameterAsText(0).replace(" ", "_")
spatialref = GetParameterAsText(1)
mouse_lon_proj = float(GetParameterAsText(2))
mouse_lat_proj = float(GetParameterAsText(3))
dem_layer = GetParameterAsText(4)
land_layer = GetParameterAsText(5)
userstream = GetParameterAsText(6) == "true"
folder_name = GetParameterAsText(7)
uploadedlayer = GetParameterAsText(8) == "true"

# Paths
directory = os.path.join(os.environ['GISHydro_DIR'])
directorygdb = os.path.join(directory, DATA_DIR, "gishydro.gdb")
demgdb = os.path.join(directory, DATA_DIR, "dem.gdb")
landusegdb = os.path.join(directory, DATA_DIR, "landuse.gdb")
soilsgdb = os.path.join(directory, DATA_DIR, "soils.gdb")

# Spatial references
sr_map = SpatialReference(SR_MAP_CODE)
sr_md = SpatialReference(SR_MD_CODE)

# Selected layers
dem_selected = os.path.join(demgdb, dem_layer)
landuse_selected = os.path.join(landusegdb, land_layer)

# Function to transform coordinates
def coord_transf(x1, y1, sr1, sr2):
    pt1 = Point(x1, y1)
    ptgeo1 = PointGeometry(pt1, sr1)
    ptgeo2 = ptgeo1.projectAs(sr2)
    pt2 = ptgeo2.lastPoint
    x2 = float(pt2.X)
    y2 = float(pt2.Y)
    return x2, y2

# Create temporary folder if it doesn't exist
temp_folder = os.path.join(directory, "projects")
os.makedirs(temp_folder, exist_ok=True)

# Create project folder
if userstream:
    project_folder = os.path.join(temp_folder, folder_name)
else:
    folder_name = time.strftime("%Y%m%d_%H%M%S") + "_" + proj_name
    project_folder = os.path.join(temp_folder, folder_name)
    os.makedirs(project_folder, exist_ok=True)

# Set environment variables
env.overwriteOutput = True
env.scratchWorkspace = project_folder
env.workspace = project_folder
env.snapRaster = os.path.join(demgdb, "neddem30")

# Get raster properties
rast = Raster(dem_selected)
cellsize = int(rast.meanCellWidth)

# Set paths based on userstream
if userstream:
    dem_selected = os.path.join(project_folder, "dem_ext.tif")
    flowdir = os.path.join(project_folder, "flowdir_ext.tif")
    infstreams = os.path.join(project_folder, "infstr_ext.tif")
else:
    hydro_gdb = os.path.join(directory, DATA_DIR, HYDRO_GDB_TEMPLATE.format(cellsize))
    flowdir = os.path.join(hydro_gdb, "flowdir{}".format(cellsize))
    infstreams = os.path.join(hydro_gdb, "infstreams{}".format(cellsize))

# Transform coordinates if necessary
if spatialref == str(SR_MAP_CODE):
    x, y = coord_transf(mouse_lon_proj, mouse_lat_proj, sr_map, sr_md)
elif spatialref == str(SR_MD_CODE):
    x = float(mouse_lon_proj)
    y = float(mouse_lat_proj)

# Create point geometry
point = PointGeometry(Point(x, y))

# Set environment extent
env.extent = Extent(x - cellsize, y - cellsize, x + cellsize, y + cellsize)

# Snap pour point
outSnapPour = SnapPourPoint(point, infstreams, 0)
pour_point_path = os.path.join(project_folder, "usr_pour_point.shp")
pour_point = RasterToPoint_conversion(outSnapPour, pour_point_path)
outSnapPour.save(os.path.join(project_folder, "pour_point.tif"))

# Get outlet coordinates
for row in SearchCursorda(pour_point_path, ["SHAPE@XY"]):
    x, y = row[0]
outletxy = management.GetCellValue(infstreams, "{} {}".format(x, y))

# Set environment extent to DEM
env.extent = Raster(dem_selected)

# Create watershed
wshed = Watershed(flowdir, os.path.join(project_folder, "pour_point.tif"))
watershed = Con(wshed >= 0, 1, IsNull(wshed))
RasterToPolygon_conversion(watershed, os.path.join(project_folder, "wshed_aux.shp"), "NO_SIMPLIFY", "VALUE")

# Set environment extent to watershed
env.extent = os.path.join(project_folder, "wshed_aux.shp")

# Dissolve and clip watershed
Dissolve_management(os.path.join(project_folder, "wshed_aux.shp"), os.path.join(project_folder, "wshed.shp"), "gridcode")
Clip_management(watershed, "", os.path.join(project_folder, "basingrid.tif"), os.path.join(project_folder, "wshed_aux.shp"), "#", "ClippingGeometry", "MAINTAIN_EXTENT")

# Project watershed and build raster attribute table
Project_management(os.path.join(project_folder, "wshed.shp"), os.path.join(project_folder, "wshed_proj.shp"), sr_map)
BuildRasterAttributeTable_management(os.path.join(project_folder, "basingrid.tif"), "Overwrite")

# Check AOI extent for land use
basingrid = Raster(os.path.join(project_folder, "basingrid.tif"))
lu_clip = Times(landuse_selected, basingrid)
lu_clip.save(os.path.join(project_folder, "landuse.tif"))
BuildRasterAttributeTable_management(os.path.join(project_folder, "landuse.tif"), "Overwrite")

# Join field
JoinField_management(os.path.join(project_folder, "landuse.tif"), "Value", landuse_selected, "VALUE", "CLASS_NAME")

# Calculate land use data percentage
try:
    lutab = SearchCursor(os.path.join(project_folder, "landuse.tif"), "", "", "Count", "")
    lucount = sum(row.getValue("Count") for row in lutab)
except:
    lucount = 0

shedtab = SearchCursor(Raster(os.path.join(project_folder, "basingrid.tif")), "", "", "Count", "")
basinarea = sum(row.getValue("Count") for row in shedtab)

ludataperc = lucount / basinarea * 100

# Check for near border watershed location
extentlyr = Clip_analysis(os.path.join(project_folder, "wshed.shp"), os.path.join(directorygdb, "boundarycheck"))
dummy_cursor = SearchCursor(extentlyr, "", "", "FID", "")
checkext = dummy_cursor.next()
extentcheck = False
if checkext is not None:
    extentcheck = True

# Delete temporary files if land use data percentage is less than 5
if ludataperc < 5:
    try:
        for filename in ["basingrid_aux.tif", "infstr.shp", "wshed_aux.shp", "infstr_aux.shp", "mask_ints.shp"]:
            Delete_management(os.path.join(project_folder, filename))
    except:
        pass

# Set parameters
SetParameterAsText(9, os.path.join(project_folder, "wshed_proj.shp"))
SetParameterAsText(10, folder_name)
SetParameterAsText(11, ludataperc)
SetParameterAsText(12, extentcheck)
