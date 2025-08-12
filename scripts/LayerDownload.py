# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 12:01:58 2023

@author: Javier.Mardones
"""

import arcpy
import os
import zipfile
import json
arcpy.env.overwriteOutput = True

# get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

projectname = arcpy.GetParameterAsText(0)
layerlist = json.loads(arcpy.GetParameterAsText(1))

optfolder = os.path.join(directory, "projects", projectname)

arcpy.env.overwriteOutput = True
arcpy.env.scratchWorkspace = optfolder
arcpy.env.workspace = optfolder

#############################
# start main proccess
#############################

layerlist_conv = {
    "wshed_dwnl": ["Watershed", "wshed"],
    "streams_dwnl": ["Streams", "nhd_strs"],
    "subshed_dwnl": ["Subwatersheds", "subshed"],
    "lfp_dwnl": ["Longest_Flow_Paths", "Longest_Path_Sub"],
    "contours_dwnl": ["Contours", "contours"],
    "dem_dwnl": ["DEM", "dem_clip"],
    "flowdir_dwnl": ["Flow_Direction", "flowdir"],
    "flowacc_dwnl": ["Flow_Accumulation", "flowacc"],
    "landuse_dwnl": ["Land_Use", "landuse"],
    "soils_dwnl": ["Soils", "soils"],
    "cn_dwnl": ["Curve_Number", "curvenumber"],
    "slope_dwnl": ["Land_Slope", "landslope"],
    "xs_dwnl": ["Cross_Sections", "xs_"]
}

zipname = projectname + "_LAYERS.zip"
# put the files that make up the shapefile into a zip file
zip = zipfile.ZipFile(os.path.join(optfolder, zipname), 'w', zipfile.ZIP_DEFLATED)

for layer in layerlist:
    for file in os.listdir(optfolder):
        if file.startswith(layerlist_conv[layer][1] + r".") or file.startswith(layerlist_conv[layer][1] + r"_"):
            try:
                zip.write(os.path.join(optfolder, file), os.path.join(layerlist_conv[layer][0], file))
            except:
                pass
zip.close()

arcpy.SetParameter(2, os.path.join(optfolder, zipname))
