# -*- coding: utf-8 -*-
"""
Modified 02/2021

@author: Javier.Mardones
"""

from arcpy import (GetParameterAsText,
                   SearchCursor,
                   UpdateCursor,
                   GetRasterProperties_management,
                   SetParameterAsText)
from arcpy.sa import Times
import os
import json
import shutil

projectname = GetParameterAsText(0)
Main_increment = float(GetParameterAsText(1))
detail = GetParameterAsText(2) == 'true'
ARC = GetParameterAsText(3)
ARF = GetParameterAsText(4) == 'true'
DelMarVa = GetParameterAsText(5)
cb_list = json.loads(GetParameterAsText(6))
prec_uSpecified = json.loads(GetParameterAsText(7))
areami2 = float(GetParameterAsText(8))
ratingtype = json.loads(GetParameterAsText(9))
minstage = json.loads(GetParameterAsText(10))
reachno = json.loads(GetParameterAsText(11))
rating = json.loads(GetParameterAsText(12))
userda = GetParameterAsText(13)
usercn = GetParameterAsText(14)
usertc = GetParameterAsText(15)
precipitation_selection = GetParameterAsText(16)

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])

optfolder = os.path.join(directory, "projects", projectname)

noaaprecipgdb = os.path.join(directory, "data","noaaprecip.gdb")
marisagdb = os.path.join(directory, "data","marisa.gdb")

# ******************************************************************************************************
# TR20 input text files with and without subwatersheds will be generated
# ******************************************************************************************************

# create "FROM NODE", "TO NODE", and "GRID CODE" lists
fn_lst = []
ai_lst = []
ln_lst = []
tn_lst = []
subriver = os.path.join(optfolder, "subrivers.shp")
subriver_len = SearchCursor(subriver, "", "", "ARCID;FROM_NODE;TO_NODE;Length", "")
for node in subriver_len:
    ai_lst.append(int(node.getValue("ARCID")))
    fn_lst.append(int(node.getValue("FROM_NODE")))
    tn_lst.append(int(node.getValue("TO_NODE")))
    ln_lst.append(int(node.getValue("Length")))
subreach_lst = [x for x in fn_lst if x in tn_lst]
reach_lst = []
len_lst = []
tolist = []
for sub in subreach_lst:
    index = fn_lst.index(sub)
    tolist.append(tn_lst[index])
    reach_lst.append(ai_lst[index])
    len_lst.append(ln_lst[index])

reach_list = []
for tn in tn_lst:
    aux = False
    for ind, fn in enumerate(fn_lst):
        if tn == fn:
            reach_list.append(ai_lst[ind])
            aux = True
    if aux == False:
        reach_list.append(0)

lst_reach = []
for rch in reach_list:
    if rch == 0:
        lst_reach.append("Outlet")
    else:
        lst_reach.append("Reach" + str(rch))

# create "slope", "area (mi^2)", "CN", and "Tc" lists using sub watersheds shapefile
sp_lst = []
da_lst = []
cn_lst = []
tc_lst = []
subshed = os.path.join(optfolder, "subshed.shp")

if userda and userda != "#":
    userdalist = json.loads(userda)
    dachange = UpdateCursor(subshed)
    count = 0
    for da in dachange:
        da.AreaMi2 = round(float(userdalist[count]),2)
        dachange.updateRow(da)
        count = count + 1

if usercn and usercn != "#":
    usercnlist = json.loads(usercn)
    cnchange = UpdateCursor(subshed)
    count = 0
    for cn in cnchange:
        cn.CurveNum = round(float(usercnlist[count]),1)
        cnchange.updateRow(cn)
        count = count + 1

if usertc and usertc != "#":
    usertclist = json.loads(usertc)
    tcchange = UpdateCursor(subshed)
    count = 0
    for t in tcchange:
        t.Tc = round(float(usertclist[count]),4)
        tcchange.updateRow(t)
        count = count + 1

ss = SearchCursor(subshed, "", "", "CurveNum;Slope;Tc;AreaMi2", "")
for att in ss:
    sp_lst.append(att.getValue("Slope"))
    da_lst.append(att.getValue("AreaMi2"))
    cn_lst.append(att.getValue("CurveNum"))
    tc_lst.append(att.getValue("Tc"))

# convert all integer/float lists into list of strings for justified text writing
ai_lst = map(str, ai_lst)
sp_lst = map(str, ["%0.2f" % x for x in sp_lst])
da_lst = map(str, ["%0.2f" % x for x in da_lst])
cn_lst = map(str, ["%0.1f" % x for x in cn_lst])
tc_lst = map(str, ["%0.2f" % x for x in tc_lst])

# write sub-area strings
sub = ""

detailstr = "NN"
if detail == True:
    detailstr = "YY"

for a, b, c, d, e in zip(ai_lst, lst_reach, da_lst, cn_lst, tc_lst):
    sub += "{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}".format("", a, b, "GAGE", c, d, e) + detailstr + "\n"

# write reach strings
NN_reach = []
for n in enumerate(len_lst):
    NN_reach.append("NN")

if len(tolist) > 0:
    outletvalue = list(set(tolist) - set(subreach_lst))[0]

toreach = []
for to,re in zip(tolist,reach_lst):
    if to == outletvalue:
        toreach.append('Outlet')
    else:
        toreach.append('Reach' + str(reach_lst[subreach_lst.index(to)]))

reachstring = ""

len_lst = ["%0.1f" % x for x in len_lst]

for h, k, l, rt in zip(toreach, reach_lst, len_lst, ratingtype):
    if rt == "XS":
        xs = "XS" + str(k)
        reachstring += "{:<10}{:<10}{:<10}{:<20}{:<20}".format("", 'Reach' + str(k), h, xs, str(l)) + "NN" + "\n"
    else:
        xs = "Struct" + str(k)
        reachstring += "{:<10}{:<10}{:<20}{:<30}".format("", 'Reach' + str(k), h, xs) + "NN" + "\n"

durlist = ["05m", "10m", "15m", "30m", "60m", "02h", "03h", "06h", "12h", "24h", "48h"]
yearlist = ["1", "2", "5", "10", "25", "50", "100", "200", "500"]

if ARF == True:
    RF_precip = []
    for y,p in zip(cb_list, prec_uSpecified):
        thecritdur = durlist[(y)%4 + 7]
        if thecritdur == "06h":
            RF_6hour  = float(p)*(1 - (0.008245*(pow(areami2,0.558))))
            RF_precip.append(RF_6hour)
        if thecritdur == "12h":
            RF_12hour = float(p)*(1 - ((0.008245/2)*(pow(areami2,0.558))) - ((0.01044/2)*(pow(areami2,0.4))))
            RF_precip.append(RF_12hour)
        if thecritdur == "24h":
            RF_24hour = float(p)*(1 - (0.01044*(pow(areami2,0.4))))
            RF_precip.append(RF_24hour)
        if thecritdur == "48h":
            RF_48hour = float(p)*(1 - (0.005*(pow(areami2,0.5169))))
            RF_precip.append(RF_48hour)
    RF_precip = ["%0.2f" %x for x in RF_precip]

else:
    RF_precip = [x for x in prec_uSpecified if x != '-']

storm_number = len(prec_uSpecified)
designstorm = []
stormid = []
stormdur = []
year = []
for i, p in zip(cb_list, prec_uSpecified):
    year.append(yearlist[i // 4])
    critdur = durlist[(i) % 4 + 7]

    stormname_map = {
        "average_opt": "p%s-%sa",
        "upper90_opt": "p%s-%su",
        "marisa41_opt": "p%s-%sm",
        "marisa42_opt": "p%s-%sm",
        "marisa81_opt": "p%s-%sm",
        "marisa82_opt": "p%s-%sm"
    }
    stormname = stormname_map.get(precipitation_selection) % (year[-1], critdur[0:2])

    stormid.append(stormname)
    stormdur.append(critdur[0:2])

    for d in durlist:
        thefilename = os.path.join(noaaprecipgdb, "p" + str(year[-1]) + "yr" + str(d) + "a")
        if precipitation_selection == "upper90_opt":
            thefilename = thefilename + "u"
        elif precipitation_selection == "marisa41_opt":
            thefilename = os.path.join(marisagdb, "r4_1_" + f"{year[-1]:03}" + "_" + d)
        elif precipitation_selection == "marisa42_opt":
            thefilename = os.path.join(marisagdb, "r4_2_" + f"{year[-1]:03}" + "_" + d)
        elif precipitation_selection == "marisa81_opt":
            thefilename = os.path.join(marisagdb, "r8_1_" + f"{year[-1]:03}" + "_" + d)
        elif precipitation_selection == "marisa82_opt":
            thefilename = os.path.join(marisagdb, "r8_2_" + f"{year[-1]:03}" + "_" + d)

        precipgrid = Times(os.path.join(optfolder, "basingrid.tif"), thefilename)        # ADDED 6/12/2020 (It doesn't make sense to not include "basingrid", otherwise its an arbitrary average)
        precavg = GetRasterProperties_management(precipgrid, "MEAN")
        precavg = float(precavg.getOutput(0))
        theavg = precavg / 1000
        designstorm.append(theavg)

###############################################
############# DESIGN STORM SCRIPT #############
###############################################

storm_response = []

tdist = [5, 10, 15, 30, 60, 120, 180, 360, 720, 1440, 2880]
tdist = [float(x)/60 for x in tdist]

stormtype = []
theprecip = []
nlines = []

imax_rd = 0
for nn in stormdur:

    iduration = int(nn)

    if iduration == 6:
        imax = 7
    elif iduration == 12:
        imax = 8
    elif iduration == 24:
        imax = 9
    elif iduration == 48:
        imax = 10

    pdist = []
    for i in range(11):
        pdist.append(float(designstorm[i + imax_rd]))

    imax_rd = imax_rd + 11

    tmax = float(iduration)
    itypeii = 0

    if pdist[0] < 0:
        itypeii = 1

    if itypeii == 0:
        stormtype.append(r'')
    else:
        stormtype.append(r'Type')

    sp = [0]*(2*(imax+1)+1)
    st = [0]*(2*(imax+1)+1)
    pdiff = []
    tdiff = []

    if itypeii == 0:
        pdiff.append(pdist[0]/2)
        tdiff.append(tdist[0]/2)
        for i in range(10):
             pdiff.append((pdist[i+1] - pdist[i])/2)
             tdiff.append((tdist[i+1] - tdist[i])/2)
        sp[imax+1] = pdist[imax]/2
        st[imax+1] = tdist[imax]/2
        for i in range(imax+1):
             sp[imax+i+2] = sp[imax + i + 1] + pdiff[i]
             st[imax+i+2] = st[imax + i + 1] + tdiff[i]
        for i in range(imax+1):
             sp[i+1] = pdist[imax] - sp[2*(imax + 1) - i]
             st[i+1] = tdist[imax] - st[2*(imax + 1) - i]

        n = int(tmax/Main_increment)

        p = []
        for i in range(n+1):
            tt = round(Main_increment * i,1)
            ii = 1
            ifound = 0
            while ifound == 0:
                ii = ii + 1
                if ii > len(st): break
                if tt >= st[ii-2] and tt < st[ii-1]:
                    ifound = ii
                    p.append(sp[ii-2] + (sp[ii-1] - sp[ii-2]) / (st[ii-1] - st[ii-2]) * (tt - st[ii-2]))

        p.append(pdist[imax])


        for i in range(n+1):
            theprecip.append(str(abs(round(p[i]/pdist[imax],4))))
        nlines.append(n+1)


storm_response.append(stormtype)
storm_response.append(theprecip)
storm_response.append(nlines)

###############################################
###############################################
###############################################

prcp = []
stormstring = ""

selection_map = {
    "average_opt": "a",
    "upper90_opt": "u",
    "marisa41_opt": "m",
    "marisa42_opt": "m",
    "marisa81_opt": "m",
    "marisa82_opt": "m"
}

noaa_str = selection_map.get(precipitation_selection, "")

for theyear, thecritdur, p, s in zip(year, stormdur, RF_precip, storm_response[0]):
    precavg = str(p)
    if s == "Type":
        stormstring = stormstring + "{:<10}{:<10}{:<20}{:<10}{:<10}".format("", "p" + theyear +"-" + thecritdur, "GAGE", precavg, "Type II") + ARC + "\n"
    else:
        stormstring = stormstring + "{:<10}{:<10}{:<20}{:<10}{:<10}".format("", "p" + theyear +"-" + thecritdur, "GAGE", precavg, "p" +theyear + "-" + thecritdur + noaa_str) + ARC + "\n"


# write rainfall distribution from designstorm text file
rainfall_string = ""
stormtype = storm_response[0]
p = storm_response[1]
n_aux = 0
for i, n in enumerate(storm_response[2]):
    rainfall_string = rainfall_string + "{:<10}{:<20}".format("", stormtype[i] + stormid[i]) + str(Main_increment) + '\n'
    for k in range(int(n/5)):
        j = [(k+1)*5-5,(k+1)*5-4,(k+1)*5-3,(k+1)*5-2,(k+1)*5-1]
        j = [x + n_aux for x in j]
        theprecip = [p[j[0]],p[j[1]],p[j[2]],p[j[3]],p[j[4]]]
        theprecip = ["%0.4f" % float(num) for num in theprecip]

        rainfall_string = rainfall_string + "{:<20}{:<10}{:<10}{:<10}{:<10}".format("", theprecip[0], theprecip[1], theprecip[2], theprecip[3]) + theprecip[4] + "\n"
    if (int(n)+1)%5 > 0:
        for k in range(int(n%5)):
            j = n - 1 + k + n_aux
            theprecip = p[j]
            theprecip = "%0.4f" % float(theprecip)
            rainfall_string = rainfall_string + "{:>26}".format(theprecip)
    rainfall_string = rainfall_string + '\n'
    n_aux = n_aux + n

# write stream cross-section block from rating table folder
rattabstring = ""
rattabstring_res = ""
for rt, ms, rn, rd in zip(ratingtype, minstage, reachno, rating):

    reach_type = rt
    reach_elev = str(round(float(ms),2))
    reach_no = str(rn)

    if reach_type == "XS":
        rattabstring = rattabstring + "{:<10}{:<10}{:<30}".format("", reach_type + reach_no, reach_elev) + "Reach" + reach_no + "\n"

        for r in rd:
            # added to format new rattab.exe output
            rattab = "{:20}{:>10}{:>10}{:>10}{:>10}{:>10}".format("", str(round(float(r[0]),2)), str(round(float(r[1]),2)), str(round(float(r[2]),2)), str(round(float(r[3]),2)), str(round(float(r[4]),2))) + "\n"
            rattabstring = rattabstring + rattab

    if reach_type == "Struct":

        rattabstring_res = rattabstring_res + "{:<10}{:<10}{:>10}".format("", reach_type + reach_no, reach_elev) + "\n"

        for r in rd:
            # added to format new rattab.exe output
            rattab = "{:<20}{:>10}{:>10}{:>10}".format("", str(round(float(r[0]),2)), str(round(float(r[1]),2)), str(round(float(r[2]),2))) + "\n"
            rattabstring_res = rattabstring_res + rattab

try:
    projname = optfolder.replace('/','\\').split('\\')[-1]
except:
    projname = ""

# Write TR20 Input file

dst = os.path.join(optfolder, "wintr20")

if os.path.exists(dst):
    for the_file in os.listdir(dst):
        file_path = os.path.join(dst, the_file)
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
            except:
                pass
else:
    os.mkdir(dst, 0o755)

inputstring = ""
inputstring = inputstring +  "{:<40}{:<10}{:<10}{:<10}".format("WinTR-20: Version 3.20", "0", "0", "1.0") + "0"
inputstring = inputstring + "\n"
inputstring = inputstring + "GISHydroWEB [" + projname + "]" + "\n"
inputstring = inputstring + "\n"
inputstring = inputstring + "SUB-AREA:" + "\n"
inputstring = inputstring + sub
inputstring = inputstring + "\n"

inputstring = inputstring + "\n"
# DelMarVa Hydrograph
if DelMarVa == "DelMarVa PRF 284":
    table_dmv = open(os.path.join(directory, "data", "lookup","Table_DMV.txt"))
    dmv_head = table_dmv.read()
    dmv_hyd = dmv_head[39:]
    inputstring = inputstring + dmv_hyd
    inputstring = inputstring + "\n"

if not reachstring == "":
    inputstring = inputstring + "STREAM REACH:" + "\n"
    inputstring = inputstring + reachstring
    inputstring = inputstring + "\n"
    inputstring = inputstring + "\n"

inputstring = inputstring + "STORM ANALYSIS:" + "\n"
inputstring = inputstring + stormstring
inputstring = inputstring + "\n"
inputstring = inputstring + "" "\n"
inputstring = inputstring + "RAINFALL DISTRIBUTION:" + "\n"
inputstring = inputstring + rainfall_string
inputstring = inputstring + "\n"

if not rattabstring == "":
    inputstring = inputstring + "STREAM CROSS SECTION:" + "\n"
    inputstring = inputstring + rattabstring
    inputstring = inputstring + "\n"

if not rattabstring_res == "":
    inputstring = inputstring + "STRUCTURE RATING:" + "\n"
    inputstring = inputstring + rattabstring_res
    inputstring = inputstring + "\n"

inputstring = inputstring + "\n"
inputstring = inputstring + "\n"
inputstring = inputstring + "GLOBAL OUTPUT:" + "\n"
inputstring = inputstring + "{:<20}{:<10}{:<10}{:<10}".format("", "1.", float(Main_increment), "YNNNN") + "YNNNNN"

defFN = os.path.join(dst, "TR20.inp")
f = open(defFN, "w")
f.write(inputstring)
f.close()

wintr20path = os.path.join(directory, "data","wintr20","WinTR20_V32.exe")

shutil.copy2(wintr20path, dst)
os.chdir(dst)    # Is this necessary?
os.system(os.path.join(dst, "WinTR20_V32.exe"))

# *******************************************************************************************************
# If TR20 error file size is greater than "0" then display a warning message to check error file
# *******************************************************************************************************
errorstring = "NA"
if os.path.exists(os.path.join(dst, "TR20.err")):

    filestats = os.path.getsize(os.path.join(dst, "TR20.err"))
    if filestats > 0:

        with open(os.path.join(dst, "TR20.err"), 'r') as file:
            errorstring = file.read()

# *******************************************************************************************************
# Change extension of "TR20in.out" file from ".out" to ".txt" and open it in text editor
# *******************************************************************************************************
outputstring = "NA"
if os.path.exists(os.path.join(dst, "TR20.out")):

    with open(os.path.join(dst, "TR20.out"), 'r') as file:
        outputstring = file.read()

if os.path.exists(os.path.join(dst, "WinTR20_V32.exe")):
    try:
        os.remove(os.path.join(dst, "WinTR20_V32.exe"))
    except:
        pass

SetParameterAsText(17, inputstring)
SetParameterAsText(18, outputstring)
SetParameterAsText(19, errorstring)
