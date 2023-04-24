import dataFlux
import numpy as np
from Globals import *
import os

def get_energy_int(f):
    f[np.where(f < 0)] = 0
    return np.sum(f, axis=1)

Nb = 101
bins = np.logspace(2, 8, Nb)
dates = ("20120907", "20130101", "20130401", "20130701", "20131001", 
            "20140101", "20140401", "20140701", "20141001",
            "20150101", "20150401", "20150701", "20151001",
            "20160101", "20160401", "20160701", "20161001",
            "20170101", "20170401", "20170701", "20171001",
            "20180101", "20180401", "20180701", "20181001",
            "20190101", "20190401", "20190701", "20191014"   )

H1 = np.zeros((Nb-1,Nb-1))
H2 = np.zeros((Nb-1,Nb-1))
for i in range(len(dates)-1):
    data = dataFlux.get_dataflux(4, dates[i], dates[i+1])
    f1 = data.get('fluxes1')
    f2 = data.get('fluxes2')
    ft = data.get('fluxesRBSP')
    f1e = get_energy_int(f1)
    f2e = get_energy_int(f2)
    fte = get_energy_int(ft)

    Ht1, xedges, yedges = np.histogram2d(f1e, fte, bins=bins)
    Ht2, xedges, yedges = np.histogram2d(f2e, fte, bins=bins)
    H1 += Ht1
    H2 += Ht2
    print("Done: "+dates[i+1])

dirname = database + "data_paper\\figure2\\"
if not os.path.isdir(dirname):
    os.makedirs(dirname)

np.savetxt(dirname+"binsX.txt",xedges)
np.savetxt(dirname+"binsY.txt",yedges)
np.savetxt(dirname+"binsZ1.txt",H1)
np.savetxt(dirname+"binsZ2.txt",H2)