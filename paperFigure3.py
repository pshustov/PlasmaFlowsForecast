import numpy as np
import datetime
import os
import dataFlux
import fitFlows

from Globals import *

def get_energy_int(f):
    f[np.where(f < 0)] = 0
    return np.sum(f, axis=1)

R = 4
date_start  = datetime.date(2012, 9, 7)
date_end    = datetime.date(2019, 10, 14)

date_current = date_start
date_step = datetime.timedelta(days=60)

SymH = np.empty(shape=(0,))
while(date_current < date_end):
    date_start_p = date_current
    date_current += date_step
    date_end_p = min(date_current, date_end)
    out_data =  dataFlux.get_dataflux(R, date_start_p.strftime('%Y%m%d'), date_end_p.strftime('%Y%m%d'), days_step=20)
    # fluxes_data =  np.concatenate( (fluxes_data, get_energy_int(out_data.get('fluxes1'))) )
    SymH = np.concatenate( (SymH, out_data.get('SymH')) )
    print('Done: %s' % (date_current))

SymH_bins = fitFlows.get_equinumeros_subests_on_param('SymH', numOfBins=40, R=4, date_start_str=date_start.strftime('%Y%m%d'), date_end_str=date_end.strftime('%Y%m%d'))


dirname = database + "data_paper\\figure3\\"
if not os.path.isdir(dirname):
    os.makedirs(dirname)

np.savetxt(dirname+'SymH.txt', SymH)
np.savetxt(dirname+'SymH_bins.txt', SymH_bins)