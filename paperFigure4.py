import numpy as np
import datetime
import os
import dataFlux
import fitFlows

from matplotlib import pyplot as plt

from Globals import *

def get_energy_bins():
    return np.array([0.0, 33.0, 54.0, 80.0,108.0,143.0,184.0,226.0,235.0,346.0,470.0,597.0,749.0,909.0,1064.0,1079.0,1575.0,1728.0,2280.0,2619.0,3618.0,4062.0,-1.0e+31,-1.0e+31,-1.0e+31,-1.0e+31])

dirname = database + "data_bins\\equals80\\"
N_fluxes = 100
N_SymH = 80
N_energy = 25

R = 4
date_start  = datetime.date(2012, 9, 7)
date_end    = datetime.date(2019, 10, 14)

date_current = date_start
date_step = datetime.timedelta(days=60)

fluxes1 = np.empty(shape=(0,25))
fluxes2 = np.empty(shape=(0,25))
SymH = np.empty(shape=(0,))
while(date_current < date_end):
    date_start_p = date_current
    date_current += date_step
    date_end_p = min(date_current, date_end)
    out_data = dataFlux.get_dataflux(R, date_start_p.strftime('%Y%m%d'), date_end_p.strftime('%Y%m%d'), days_step=20)
    fluxes1 = np.concatenate( (fluxes1, out_data.get('fluxes1')) )
    fluxes2 = np.concatenate( (fluxes2, out_data.get('fluxes2')) )
    SymH = np.concatenate( (SymH, out_data.get('SymH')) )
    print('Done: %s' % (date_current))

# SymH_bins = fitFlows.get_equinumeros_subests_on_param('SymH', numOfBins=N_SymH, R=4, date_start_str=date_start.strftime('%Y%m%d'), date_end_str=date_end.strftime('%Y%m%d'))
SymH_bins = np.linspace(-200, 50, N_SymH+1)

indsH = np.digitize(SymH,SymH_bins)
fluxes_bins = np.logspace(-3,8,N_fluxes+1)

flux_hist = np.zeros((N_fluxes, N_energy, N_SymH))
for i in range(N_SymH):
    print(i/N_SymH)
    for j in range(N_energy):
        flux_hist[:,j,i], bins = np.histogram(fluxes1[indsH==(i+1),j], bins=fluxes_bins)


if not os.path.isdir(dirname):
    os.makedirs(dirname)

np.savetxt(dirname+'energy_bins.txt', get_energy_bins())   
np.savetxt(dirname+'fluxes_bins.txt', fluxes_bins)   
np.savetxt(dirname+'SymH_bins.txt', SymH_bins)    


dirname = dirname + "flux_hist\\"
if not os.path.isdir(dirname):
    os.makedirs(dirname)

for i in range(N_energy):
    np.savetxt(dirname+'flux_hist_'+str(i)+'.txt', flux_hist[:,i,:])
