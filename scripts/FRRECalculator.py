# -*- coding: utf-8 -*-
"""
Modified 02/2021

@author: Javier.Mardones
"""

#  ADJUSTED SCRIPT FROM GISHYDRONXT TO GISHYDROWEB
#  By: JAVER MARDONES
#  Last Modified: 05/2022
#
#  ---------------------------------------------------------------------------
#
#  Original description:
#
#  ---------------------------------------------------------------------------
#
# Program to estimate flood frequency in Maryland
#

#
#  This program calculates the fixed region Flood Frequency based on the
#  regression equations developed by Wilbert Thomas as part of the study
#  conducted by Moglen, Thomas, and Cuneo (2006).
#
#  Program originally developed by Gary Tasker in 1999 to accompany the
#  regression equations developed by Dillow (1996).
#
#  Program modified by Glenn Moglen in 2010, and 2016 to accompany the fixed region
#  regression equations developed by Thomas and Moglen (2010, 2015).
#

# Import Libraries

from arcpy import (GetParameterAsText, SetParameterAsText)
import numpy as np
import os
import json

####
global set_i
set_i = 0
region = GetParameterAsText(set_i); set_i = set_i + 1
area = float(GetParameterAsText(set_i)); set_i = set_i + 1
parameter1 = GetParameterAsText(set_i); set_i = set_i + 1
parameter2 = GetParameterAsText(set_i); set_i = set_i + 1
gageid = GetParameterAsText(set_i); set_i = set_i + 1
eqyr = GetParameterAsText(set_i); set_i = set_i + 1


if eqyr == "2020":
    gagefile = os.path.join("2020", "gagerat2020.txt")
    regfile = os.path.join("2020", region + "_frredat20.txt")
elif eqyr == "2022":
    gagefile = os.path.join("2022", "gagerat2022.txt")
    regfile = os.path.join("2022", region + "_frredat22.txt")

#get path from environment variable
directory = r"" + str(os.environ['GISHydro_DIR'])
global directorytasker
directorytasker = os.path.join(directory, "data", "frredat")

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


### REGRESSION EQUATION CALCULATOR

def REG(filedir, paramvals, gageid):
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

    paramvals = [float(i) for i in paramvals[:npred-1]]
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


    global set_i
    SetParameterAsText(set_i, json.dumps(estim_par)); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps(warning_message)); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps([[round(j,0) for j in i] for i in cl])); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps([[round(j,0) for j in i] for i in cu])); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps([round(i,0) for i in yhat_list])); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps([round(i,1) for i in sepc_list])); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps([round(i,2) for i in eqyrs_list])); set_i = set_i + 1
    SetParameterAsText(set_i, json.dumps([round(i,4) for i in sepred_list])); set_i = set_i + 1

### Check gage list

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

filedir = os.path.join(directorytasker, regfile)
paramvals = [area, parameter1, parameter2]

REG(filedir, paramvals, gageid)