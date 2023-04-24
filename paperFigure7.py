import numpy as np
import datetime
import os
import dataFlux
import fitFlows

from matplotlib import pyplot as plt

from Globals import *

def get_energy_bins():
    return np.array([0.0, 33.0, 54.0, 80.0,108.0,143.0,184.0,226.0,235.0,346.0,470.0,597.0,749.0,909.0,1064.0,1079.0,1575.0,1728.0,2280.0,2619.0,3618.0,4062.0,-1.0e+31,-1.0e+31,-1.0e+31,-1.0e+31])

dirname = database + "data_paper\\figure7\\"
N_fluxes = 100
N_SymH = 80
N_energy = 25

R = 4
date_start  = datetime.date(2016, 2, 12)
date_end    = datetime.date(2016, 11, 14)

date_current = date_start
date_step = datetime.timedelta(days=60)

times = np.empty(shape=(0,))
fluxes1 = np.empty(shape=(0,25))
fluxes2 = np.empty(shape=(0,25))
SymH = np.empty(shape=(0,))
while(date_current < date_end):
    date_start_p = date_current
    date_current += date_step
    date_end_p = min(date_current, date_end)
    out_data = dataFlux.get_dataflux(R, date_start_p.strftime('%Y%m%d'), date_end_p.strftime('%Y%m%d'), days_step=20)
    times = np.concatenate( (times, out_data.get('times')) )
    fluxes1 = np.concatenate( (fluxes1, out_data.get('fluxes1')) )
    fluxes2 = np.concatenate( (fluxes2, out_data.get('fluxes2')) )
    SymH = np.concatenate( (SymH, out_data.get('SymH')) )
    print('Done: %s' % (date_current))

if not os.path.isdir(dirname):
    os.makedirs(dirname)

# SymH = SymH[::10]
# fluxes1 = fluxes1[::10,:]

np.savetxt(dirname+'times.txt', times)   
np.savetxt(dirname+'SymH.txt', SymH)   
np.savetxt(dirname+'fluxes.txt', fluxes1, fmt='%1.4e')
np.savetxt(dirname+'energy_bins.txt', get_energy_bins())   