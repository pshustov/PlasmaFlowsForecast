import numpy as np
import os

from matplotlib import pyplot as plt

from Globals import *

dirname = database + "data_paper\\figure4\\"
energy_bins = np.loadtxt(dirname+'energy_bins.txt')
fluxes_bins = np.loadtxt(dirname+'fluxes_bins.txt')
SymH_bins   = np.loadtxt(dirname+'SymH_bins.txt')


N_energy = len(energy_bins) - 1
N_fluxes = len(fluxes_bins) - 1
N_SymH   = len(SymH_bins) - 1


flux_hist = np.zeros((N_fluxes, N_SymH, N_energy))
dirname = database + "data_paper\\figure4\\flux_hist\\"
for i in range(N_energy):
    flux_hist[:,:,i] = np.loadtxt(dirname+'flux_hist_'+str(i)+'.txt')


fig = plt.figure(figsize=[1.5*6.5,1.5*19.5])
fig_dirname = database+'figures\\step_fit\\test1\\'
if not os.path.isdir(fig_dirname):
    os.makedirs(fig_dirname)
paramBinsNames = [('%.0f..%.0f'%(par1,par2)) for par1,par2 in zip(SymH_bins[:-1],SymH_bins[1:])]
paramBinsNames.append('SymH')


for i in range(N_energy):
    binsN = flux_hist[:,:,i]
    norm = np.sum(binsN,axis=0)
    binsData = (binsN / norm).T

    fig_name = ('figure_distrib_fit_E=%gkeV.png' % energy_bins[i+1])
    ax = fig.add_axes([0.13, 0.02, 0.85, 0.94])
    
    for j in range(binsData.shape[0]):
        norm_binsdata = np.max(binsData[j,:]) / 0.95
        norm2 = 3
        stepY = np.log10(binsData[j,:]/norm_binsdata) / norm2 + j + 1
        ax.step(fluxes_bins[1:], stepY, where='pre', color='tab:blue')
        ax.axhline(j, color='gray', linewidth=0.5)
        ax.set_yticks(np.arange(len(paramBinsNames)))
        ax.set_yticklabels(paramBinsNames)
        ax.set_xscale('log')
    
    ax.set_title('%g keV' % energy_bins[i+1])
    ax.set_ylim([0,binsData.shape[0]+0.5])
    ax.set_xlim([np.min(fluxes_bins), np.max(fluxes_bins)])
    fig.savefig(fig_dirname + fig_name)
    fig.clf()

plt.close(fig)
