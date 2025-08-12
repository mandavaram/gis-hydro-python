# -*- coding: utf-8 -*-
"""
Modified 05/2024

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   env,
                   GetRasterProperties_management,
                   SetParameterAsText)
from arcpy.sa import Times
import json
import os

projectname = GetParameterAsText(0)
cb_list = json.loads(GetParameterAsText(1))
precipitation_selection  = GetParameterAsText(2)

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

optfolder = os.path.join(directory, "projects", projectname)

noaaatlasgdb = os.path.join(directory, "data","noaaprecip.gdb")
marisagdb = os.path.join(directory, "data","marisa.gdb")

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = os.path.join(optfolder, "dem_clip.tif")

# ******************************************************************************************************
# Use selected duration and year to compute average precipitation
# ******************************************************************************************************

yearlist = ["1", "2", "5", "10", "25", "50", "100", "200", "500"]
durlist = ["06", "12", "24", "48"]

thecritavg = []
year = []
duration = []
for cb in cb_list:

    # following year and duration list is to create avg prec list
    theyear = int(yearlist[cb // 4])  # "//" will floor the value to get respective indexed year from above
    thecritdur = durlist[cb % 4]
    year.append(theyear)
    duration.append(thecritdur)

    file_map = {
        "upper90_opt": os.path.join(noaaatlasgdb, "p" + str(theyear) + "yr" + thecritdur + "hau"),
        "marisa41_opt": os.path.join(marisagdb, "r4_1_" + f"{theyear:03}" + "_" + thecritdur + "h"),
        "marisa42_opt": os.path.join(marisagdb, "r4_2_" + f"{theyear:03}" + "_" + thecritdur + "h"),
        "marisa81_opt": os.path.join(marisagdb, "r8_1_" + f"{theyear:03}" + "_" + thecritdur + "h"),
        "marisa82_opt": os.path.join(marisagdb, "r8_2_" + f"{theyear:03}" + "_" + thecritdur + "h")
    }
    
    thefilename = file_map.get(precipitation_selection, os.path.join(noaaatlasgdb, "p" + str(theyear) + "yr" + thecritdur + "ha"))

    basingrid = Times(os.path.join(optfolder, "basingrid.tif"), thefilename)
    precavg = GetRasterProperties_management(basingrid, "MEAN")
    precavg = float(precavg.getOutput(0))

    thecritavg.append(round(precavg/1000, 2))


SetParameterAsText(3, json.dumps(thecritavg));
SetParameterAsText(4, json.dumps(year))
SetParameterAsText(5, json.dumps(duration))
