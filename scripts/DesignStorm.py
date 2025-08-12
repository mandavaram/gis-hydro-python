"""
Created by Glen Moglen and edited by Javier Mardones

Last modification Modified 05/2024

"""

from arcpy import GetParameterAsText, SetParameterAsText
import json

stormdur = json.loads(GetParameterAsText(0))
stormdata = json.loads(GetParameterAsText(1))
dt = float(GetParameterAsText(2))

tdist = [x/60 for x in [5, 10, 15, 30, 60, 120, 180, 360, 720, 1440, 2880]]

stormtype = []
theprecip = []
nlines = []

imax_rd = 0
for nn in stormdur:
    iduration = int(nn)
    imax = {6: 7, 12: 8, 24: 9, 48: 10}[iduration]

    pdist = [float(stormdata[i + imax_rd]) for i in range(11)]
    imax_rd += 11

    tmax = float(iduration)
    itypeii = int(pdist[0] < 0)

    stormtype.append(r'rt' if itypeii == 0 else r'Type')

    sp = [0]*(2*(imax+1)+1)
    st = [0]*(2*(imax+1)+1)
    pdiff = [pdist[0]/2] + [(pdist[i+1] - pdist[i])/2 for i in range(10)]
    tdiff = [tdist[0]/2] + [(tdist[i+1] - tdist[i])/2 for i in range(10)]

    sp[imax+1] = pdist[imax]/2
    st[imax+1] = tdist[imax]/2
    for i in range(imax+1):
        sp[imax+i+2] = sp[imax + i + 1] + pdiff[i]
        st[imax+i+2] = st[imax + i + 1] + tdiff[i]
        sp[i+1] = pdist[imax] - sp[2*(imax + 1) - i]
        st[i+1] = tdist[imax] - st[2*(imax + 1) - i]

    n = int(tmax/dt)

    p = [sp[ii-2] + (sp[ii-1] - sp[ii-2]) / (st[ii-1] - st[ii-2]) * (round(dt * i,1) - st[ii-2]) 
         for i in range(n+1) 
         for ii in range(1, len(st)) 
         if round(dt * i,1) >= st[ii-2] and round(dt * i,1) < st[ii-1]]

    p.append(pdist[imax])

    theprecip.extend([str(abs(round(p[i]/pdist[imax],4))) for i in range(n+1)])
    nlines.append(n+1)

SetParameterAsText(3, json.dumps(stormtype))
SetParameterAsText(4, json.dumps(theprecip))
SetParameterAsText(5, json.dumps(nlines))
