
# -*- coding: utf-8 -*-
"""
Modified 08/2021

@author: Javier.Mardones
"""
from arcpy import (env,
                   GetParameterAsText,
                   SearchCursor,
                   Raster,
                   SetParameterAsText)
from arcpy.sa import ZonalStatisticsAsTable
import os
import json
import numpy as np

### Input

set_i = 0
projectname = GetParameterAsText(set_i); set_i = set_i + 1
gageid = GetParameterAsText(set_i); set_i = set_i + 1
landslope_user = GetParameterAsText(set_i); set_i = set_i + 1
imp_user = GetParameterAsText(set_i); set_i = set_i + 1
asoil_user = GetParameterAsText(set_i); set_i = set_i + 1
lime_user = GetParameterAsText(set_i); set_i = set_i + 1
reg_version = GetParameterAsText(set_i); set_i = set_i + 1

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])
directorygdb = os.path.join(directory, "data","gishydro.gdb")
global directorytasker
directorytasker = os.path.join(directory,"data","frredat")

optfolder = os.path.join(directory, "projects", projectname)

# paths
basingrid = os.path.join(optfolder, "basingrid.tif")
prov = os.path.join(directorygdb, "md_phyisio_provinces")

# Environmental variables

env.overwriteOutput = True
env.scratchWorkspace = optfolder
env.workspace = optfolder
env.snapRaster = basingrid

# *******************************************************************************************************
# don"t add "theVTab" to TOC -- it will change list by drawing order to list
# by source which will prohibit addition of new layers
# *******************************************************************************************************
ZonalStatisticsAsTable(prov, "PROVINCE", basingrid, "theVTab", "DATA", "ALL")

# *******************************************************************************************************
sumarea = 0
theVTab = SearchCursor("theVTab", "", "", "Count", "")
for each in theVTab:
    count = each.getValue("Count")
    sumarea = sumarea + count
sumArea = sumarea
del each

# ******************************************************************************************************
# define initial values and lists. add and index those lists as part of dictionary
# ******************************************************************************************************

Q1p25 = 0
Q1p50 = 0
Q2 = 0
Q5 = 0
Q10 = 0
Q25 = 0
Q50 = 0
Q100 = 0
Q200 = 0
Q500 = 0
Q1p25list = [0, 0, 0, 0, 0, 0, 0, 0]
Q1p50list = [0, 0, 0, 0, 0, 0, 0, 0]
Q2list = [0, 0, 0, 0, 0, 0, 0, 0]
Q5list = [0, 0, 0, 0, 0, 0, 0, 0]
Q10list = [0, 0, 0, 0, 0, 0, 0, 0]
Q25list = [0, 0, 0, 0, 0, 0, 0, 0]
Q50list = [0, 0, 0, 0, 0, 0, 0, 0]
Q100list = [0, 0, 0, 0, 0, 0, 0, 0]
Q200list = [0, 0, 0, 0, 0, 0, 0, 0]
Q500list = [0, 0, 0, 0, 0, 0, 0, 0]

qlist = {"1": [], "2": [], "3": [], "4": [], "5": [], "6": [], "7": [], "8": [], "9": [], "10": []}
qlist["1"].extend(Q1p25list)
qlist["2"].extend(Q1p50list)
qlist["3"].extend(Q2list)
qlist["4"].extend(Q5list)
qlist["5"].extend(Q10list)
qlist["6"].extend(Q25list)
qlist["7"].extend(Q50list)
qlist["8"].extend(Q100list)
qlist["9"].extend(Q200list)
qlist["10"].extend(Q500list)

out_rast = Raster(basingrid)
cellsize = out_rast.meanCellWidth
cellsq = cellsize * cellsize

shedtab = SearchCursor(basingrid, "", "", "Count", "")
for row in shedtab:
    basinarea = row.getValue("Count")
areami2 = float((basinarea * cellsq) / 2588881)  # conversion into sq miles

discharge_values = []
flood_intervals = []
tasker_all = []
regioncount = 0

# loop to begin tasker handling and area weighted analysis
theVTab = SearchCursor("theVTab", "", "","Province;Count", "")
for row in theVTab:
    AreaField = float(row.getValue("Count"))
    areapercent = float((AreaField / sumArea) * 100)

    intasker_sub = []
    if row.getValue("Province") == "A":
        region = 'A'
    elif row.getValue("Province") == "B" or row.getValue("Province") == "P":
        region = 'P'
    elif row.getValue("Province") == "W":
        region = 'WC'
    elif row.getValue("Province") == "E":
        region = 'EC'

    tasker_response = []

    def MLTPLY(X,Y,K1,K2):

    # --------------------------------------------------------------
    #  X IS K1 VECTOR
    #  Y IS K1*K2 MATRIX
    #  PROD = X*Y IS A K1*K2 MATRIX
    # --------------------------------------------------------------
        PROD = []
        for k in range(K2):
            var_sum = 0
            for j in range(K1):
                try:
                    var_sum = var_sum + X[j]*Y[k][j]
                except:
                    var_sum = var_sum + X[j]*Y[j]
            PROD.append(var_sum)
        return PROD


    def EYEARS(sig,zp,skew,vpi):

        #  compute equivalent years
        #
        #  Wilson-Hilferty approximation (Handbook of Hydrology, eqn 18.2.29)
        whkp=(2./skew)*(1.+(skew*zp/6.0)-skew**2/36.)**3-2.0/skew
        #  Handbook of hydrology eqn 18.4.11
        hard=1.0+skew*whkp+0.5*(1.+.75*skew**2)*whkp**2
        eqyrs = sig**2*(hard)/vpi
        return eqyrs


    def GADJ(area, gageid, eqyrs, ip, region):
    #
    #  subroutine adjusts estimates for nearby gaging station
    #  ref. Sauer, 1974, USGS WRI 52-73.
    #

        rat_txt = open(os.path.join(directorytasker, gagefile), "r")

        # Splits the element by "\n"
        rat_lines = rat_txt.readlines()
        rat_txt.close()
        rat_lines = [line[:-1] for line in rat_lines]

        gagepr = []
        sid = []
        agage = []
        r = []
        e = []
        for line in rat_lines:
            aux_line = line.split()
            gagepr.append(aux_line[0])
            sid.append(aux_line[1])
            agage.append(aux_line[2])
            r.append(aux_line[3::2])
            e.append(aux_line[4::2])

        radj=1.0
        yadj=0.0
        iflag_aux = 0
        for i in range(len(rat_lines)):
            if sid[i] == str(gageid) and gagepr[i].lower() == region.lower():
                iflag_aux = 1
                delta = abs(area-float(agage[i]))
                iflag = 3
                if delta < .5*float(agage[i]):
                    radj=float(r[i][ip])-delta*(float(r[i][ip])-1.)/(0.5*float(agage[i]))
                    yadj=((float(e[i][ip])-eqyrs)/(float(r[i][ip])-1.0))*(radj-1.0)
                else:
                    iflag = 1
            elif iflag_aux == 0:
                iflag=2

        return radj,yadj,iflag


    def REG(filedir, paramvals, gageid, areapercent):
        mrd_txt = open(filedir, "r")

        # Splits the element by "\n"
        mrd_lines = mrd_txt.readlines()
        mrd_txt.close()

        line0 = mrd_lines[0]
        region_name = str(line0.strip())

        line1 = mrd_lines[1].split()
        region = str(line1[0])
        npred = int(line1[1])
        skew = float(line1[2])

        line2 = mrd_lines[2].split()
        stut = [float(i) for i in line2]

        line3 = mrd_lines[3].split()
        ak = [float(i) for i in line3]

        line4 = mrd_lines[4].split()
        vmodel = [float(i)**2 for i in line4]

        line5 = mrd_lines[5].split(",")
        parameters = [str(i.strip()) for i in line5]

        line6 = mrd_lines[6].split(',')
        partype = [str(i.strip()) for i in line6]

        line7 = mrd_lines[7].split()
        minrange = [float(line7[i]) for i in range(len(line7)) if i % 2 == 0]
        maxrange = [float(line7[i]) for i in range(len(line7)) if i % 2 == 1]

        bs = []
        xt = []
        for j in range(10):
            bs.append([float(i) for i in mrd_lines[8 + j].split()])

            xt_aux = []
            for k in range(npred):
                xt_aux.append([float(i) for i in mrd_lines[18 + k + npred*j].split()])
            xt.append(xt_aux)

        iflag = 0
        if not gageid:
            gage_name = "No Adjustment"
        else:
            gage_name = gageid
            iflag = 3


        params = ["Region"] + parameters + ["Skew", "Gage ID"]
        params2 = [region_name] + paramvals + [skew, gage_name]
        estim_par = [params, params2]

        v = [1.0]
        for i,ptype in zip(paramvals,partype):
            if ptype == "log":
                v.append(np.log10(float(i)))
            elif ptype == "log+1":
                v.append(np.log10(float(i+1)))
            elif ptype == "nolog":
                v.append(float(i))

        vt = v
        ivpi = 0
        cu = []
        cl = []
        yhat_list = []
        sepc_list = []
        eqyrs_list = []
        sepred_list = []

        for ip in range(10):
            yhat = sum([x*y for x, y in zip(bs[ip], v)])
            yhat = 10**yhat

            # Compute CI

            xtxi = []
            for i in range(npred):
                xtxi_aux = []
                for j in range(npred):
                    xtxi_aux.append(xt[ip][j][i])
                xtxi.append(xtxi_aux)

            temp = MLTPLY(v,xtxi,npred,npred)
            temp2 = MLTPLY(temp,vt,npred,1)
            temp2 = float(temp2[0])

            if temp2 < 0:
                temp2 = 0
                ivpi = ivpi + 1

            vpi = vmodel[ip] + temp2
            sepred = np.sqrt(vpi)

            # compute equivalent years

            # Regional estimate of sigma
            sig = 0.2353
            eqyrs = EYEARS(sig,ak[ip],skew,vpi)

            if iflag == 3:
                (radj,yadj,iflag) = GADJ(paramvals[0], gageid, eqyrs, ip, region)
                yhat = 10**(np.log10(yhat)) * radj

                vpi = vpi*eqyrs/(eqyrs+yadj)
                sepred = np.sqrt(vpi)
                eqyrs = eqyrs+yadj

            sepc = 100*np.sqrt(np.exp(vpi*5.302)-1)

            cu_aux = []
            cl_aux = []
            for iis in range(4):
                t = 10**(stut[iis]*vpi**.5)
                cu_aux.append(yhat*t)
                cl_aux.append(yhat/t)
            cu.append(cu_aux)
            cl.append(cl_aux)

            yhat_list.append(yhat)
            sepc_list.append(sepc)
            eqyrs_list.append(eqyrs)
            sepred_list.append(sepred)

        warning_message = []
        if iflag == 3:
            if gageid == "01591500":
                warning_message.append("Estimates adjusted for proximity to stations 01591400 and 01591500")
            elif gageid == "01588500":
                warning_message.append("Estimates adjusted for proximity to stations 01588500 and 01589000")
            else:
                warning_message.append("Estimates adjusted for proximity to station %s" % (gageid))

        if iflag == 2:
            warning_message.append("Station, %s, unknown: NO ADJUSTMENT MADE" % (gageid))

        if not gageid:
            warning_message.append("No station selected: NO ADJUSTMENT MADE")

        if iflag == 1:
            warning_message.append("Difference in drainage area for Station %s too great: NO ADUSTMENT MADE" % (gageid))

        if ivpi > 0:
            warning_message.append("Warning: VPI is negative %s" % (ivpi))

        for i, j, pvalue, pname in zip(minrange, maxrange, paramvals, parameters):

            if pvalue < i or pvalue > j:
                warning_message.append("WARNING - " + pname + " out of range of observed data")

        tasker_response.append(estim_par)
        tasker_response.append(warning_message)
        tasker_response.append([[round(j*areapercent/100,0) for j in i] for i in cl])
        tasker_response.append([[round(j*areapercent/100,0) for j in i] for i in cu])
        tasker_response.append([round(i*areapercent/100,0) for i in yhat_list])
        tasker_response.append([round(i,1) for i in sepc_list])
        tasker_response.append([round(i,2) for i in eqyrs_list])
        tasker_response.append([round(i,4) for i in sepred_list])

    ### Check gage list

    if reg_version == "2020":
        gagefile = os.path.join("2020", "gagerat2020.txt")
        regfile = os.path.join("2020", region + "_frredat20.txt")
    elif reg_version == "2022":
        gagefile = os.path.join("2022", "gagerat2022.txt")
        regfile = os.path.join("2022", region + "_frredat22.txt")

    rat_txt = open(os.path.join(directorytasker, gagefile), "r")
    rat_lines = rat_txt.readlines()
    rat_txt.close()
    rat_lines = [line[:-1] for line in rat_lines]

    gagelist = []
    for line in rat_lines:
        aux_line = line.split()
        gagelist.append(aux_line[1])

    if not gageid in gagelist:
        gageid = False

    ### Run FRRE Computations

    if region == "A":
        filedir = os.path.join(directorytasker, regfile)
        paramvals = [round(float(areami2),3), round(float(landslope_user),5)]
    elif region == "P":
        filedir = os.path.join(directorytasker, regfile)
        paramvals = [round(float(areami2),3), round(float(imp_user),2),round(float(lime_user),2)]
    elif region == "WC":
        filedir = os.path.join(directorytasker, regfile)
        paramvals = [round(float(areami2),3), round(float(imp_user),2), round(float(asoil_user),2)]
    elif region == "EC":
        filedir = os.path.join(directorytasker, regfile)
        paramvals = [round(float(areami2),3), round(float(landslope_user)*100,3), round(float(asoil_user),2)]

    REG(filedir, paramvals, gageid, areapercent)

    regioncount = regioncount + 1
    tasker_all.append(tasker_response)

    ##############################
    ########Area Weigth###########
    ##############################

    yhat_Q = tasker_response[4]

    Q1p25 = Q1p25 + float(yhat_Q[0])
    Q1p50 = Q1p50 + float(yhat_Q[1])
    Q2 = Q2 + float(yhat_Q[2])
    Q5 = Q5 + float(yhat_Q[3])
    Q10 = Q10 + float(yhat_Q[4])
    Q25 = Q25 + float(yhat_Q[5])
    Q50 = Q50 + float(yhat_Q[6])
    Q100 = Q100 + float(yhat_Q[7])
    Q200 = Q200 + float(yhat_Q[8])
    Q500 = Q500 + float(yhat_Q[9])

    # *****************************************************************************
    # compute discharge and assign confidence intervals to qlist entries
    # *****************************************************************************

    cl = tasker_response[2]
    cu = tasker_response[3]
    for i in range(10):
        c1 = [round((float(cl[i][0]) * areapercent) / 100,0)]
        c2 = [round((float(cu[i][0]) * areapercent) / 100,0)]
        c3 = [round((float(cl[i][1]) * areapercent) / 100,0)]
        c4 = [round((float(cu[i][1]) * areapercent) / 100,0)]
        c5 = [round((float(cl[i][2]) * areapercent) / 100,0)]
        c6 = [round((float(cu[i][2]) * areapercent) / 100,0)]
        c7 = [round((float(cl[i][3]) * areapercent) / 100,0)]
        c8 = [round((float(cu[i][3]) * areapercent) / 100,0)]
        qlist = [c1, c2, c3, c4, c5, c6, c7, c8]
        flood_intervals.append(qlist)


lists = {i: [el[0] for el in v] for i, v in enumerate(flood_intervals, start=1)}

# ******************************************************************************************************
# index and count number of sub-lists -- prepare for TaskerString function
# ******************************************************************************************************

#*******************************************************************************************************
# Count number of sub-lists in a list -- special case for multiple provinces
#*******************************************************************************************************
def FloodIntervalLists(lists):
        count = 0
        for item in lists:
            if isinstance(item,list):
                count += 1 + FloodIntervalLists(item)
            else:
                count += 1
        return count
number = FloodIntervalLists(lists)

if number > 10:
    q_list1 = [sum(i) for i in zip(lists[1], lists[11])]
    q_list2 = [sum(i) for i in zip(lists[2], lists[12])]
    q_list3 = [sum(i) for i in zip(lists[3], lists[13])]
    q_list4 = [sum(i) for i in zip(lists[4], lists[14])]
    q_list5 = [sum(i) for i in zip(lists[5], lists[15])]
    q_list6 = [sum(i) for i in zip(lists[6], lists[16])]
    q_list7 = [sum(i) for i in zip(lists[7], lists[17])]
    q_list8 = [sum(i) for i in zip(lists[8], lists[18])]
    q_list9 = [sum(i) for i in zip(lists[9], lists[19])]
    q_list10 = [sum(i) for i in zip(lists[10], lists[20])]
else:
    q_list1 = lists[1]
    q_list2 = lists[2]
    q_list3 = lists[3]
    q_list4 = lists[4]
    q_list5 = lists[5]
    q_list6 = lists[6]
    q_list7 = lists[7]
    q_list8 = lists[8]
    q_list9 = lists[9]
    q_list10 = lists[10]
q_list_all = [q_list1, q_list2, q_list3, q_list4, q_list5, q_list6, q_list7, q_list8, q_list9, q_list10]
it_values = ['1.25', '1.50', '2', '5', '10', '25', '50', '100', '200', '500']

# ******************************************************************************************************
# discharge computation based on province
# ******************************************************************************************************

qcfs = [Q1p25,Q1p50,Q2,Q5,Q10,Q25,Q50,Q100,Q200,Q500]

#############
## OUTPUTS ##
#############

SetParameterAsText(set_i, json.dumps(it_values)); set_i = set_i + 1
SetParameterAsText(set_i, json.dumps(q_list_all)); set_i = set_i + 1
SetParameterAsText(set_i, json.dumps(qcfs)); set_i = set_i + 1
SetParameterAsText(set_i, regioncount); set_i = set_i + 1
SetParameterAsText(set_i, json.dumps(tasker_all)); set_i = set_i + 1


