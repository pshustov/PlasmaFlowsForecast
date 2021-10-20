import pickle
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyspedas
import os


from scipy.optimize.minpack import curve_fit
from projectionAtR import load_projectionAt
from projectionAtR import load_binsN
from projectionAtR import fit_fluxes


def plot2dHyst(fig,ax, xBins, yBins, zData):
    X,Y = np.meshgrid(xBins,yBins)
    ph = ax.pcolor(X,Y,zData,shading='flat',norm=matplotlib.colors.LogNorm())
    fig.colorbar(ph, ax=ax)

def makeplot_fluxes_distrib_onDST(fig, ax, date_start, date_end, R, yName, zName, fBin, Nx=40, Ny=40):
    xBins = np.logspace(0, 6, Nx+1)
    yBins = np.linspace(-200,50,Ny+1)

    N_load_step, N_load_step_part = 60, 20
    binsN = load_binsN(date_start, date_end, R, xBins, yBins, fBin, yName, zName, N_load_step, N_load_step_part)
    binsData = (binsN.T / np.sum(binsN,axis=1)).T

    plot2dHyst(fig,ax, xBins, yBins, binsData)
    ax.set_xscale('log')
    
def fit_fluxes_onesE(fluxBins, binsN):
    binsData = (binsN.T / np.sum(binsN,axis=1)).T
    xs = np.log10(fluxBins[1:])
    dx = np.mean(np.diff(xs))
    f = lambda x, X, s: np.exp(-(x-X)**2 / (2*s*s)) / np.sqrt(2*np.pi*s*s) * dx
    popt = np.zeros((binsData.shape[0],2))

    for i in range(binsData.shape[0]):
        p0 = [3, 1]
        popt[i,:], pcov = curve_fit(f, xs, binsData[i,:])

    return popt


def makefigure(date_start, date_end, R, Nx, Ny):
    fig, axs = plt.subplots(3,2)
    makeplot_fluxes_distrib_onDST(fig, axs[0,0], date_start, date_end, R, yName='SymH', zName='fluxes1', fBin=[3], Nx=Nx, Ny=Ny)
    makeplot_fluxes_distrib_onDST(fig, axs[0,1], date_start, date_end, R, yName='SymH', zName='fluxes2', fBin=[3], Nx=Nx, Ny=Ny)
    makeplot_fluxes_distrib_onDST(fig, axs[1,0], date_start, date_end, R, yName='SymH', zName='fluxes1', fBin=[9,10], Nx=Nx, Ny=Ny)
    makeplot_fluxes_distrib_onDST(fig, axs[1,1], date_start, date_end, R, yName='SymH', zName='fluxes2', fBin=[9,10], Nx=Nx, Ny=Ny)
    makeplot_fluxes_distrib_onDST(fig, axs[2,0], date_start, date_end, R, yName='SymH', zName='fluxes1', fBin=[12,14], Nx=Nx, Ny=Ny)
    makeplot_fluxes_distrib_onDST(fig, axs[2,1], date_start, date_end, R, yName='SymH', zName='fluxes2', fBin=[12,14], Nx=Nx, Ny=Ny)

    axs[0,0].set_title('North fluxes')
    axs[0,0].set_ylabel('E = 100keV\nSymH')
    axs[1,0].set_ylabel('E = 500keV\nSymH')
    axs[2,0].set_ylabel('E = 1MeV\nSymH')
    axs[2,0].set_xlabel('Fluxes')
    axs[0,0].set_xticklabels([])
    axs[1,0].set_xticklabels([])

    axs[0,1].set_title('South fluxes')
    axs[2,1].set_xlabel('Fluxes')
    axs[0,1].set_xticklabels([])
    axs[0,1].set_yticklabels([])
    axs[1,1].set_xticklabels([])
    axs[1,1].set_yticklabels([])
    axs[2,1].set_yticklabels([])
    plt.subplots_adjust(left=0.09, bottom=0.05, right=0.98, top=0.95, wspace=0, hspace=0.1)
    plt.show()

def get_energy(isNotNaN=False):
    mageis_vars = pyspedas.rbsp.mageis(['2019-01-01/00:00:00', '2019-01-01/01:00:00'], probe='a', level='l3', rel='rel04', notplot=True)
    energy = mageis_vars.get('FEDU').get('v2') 
    energy = energy[np.where(energy>0)]
    return energy


def matlabfun():
    date_start = datetime.date(2012, 9, 7)
    date_end = datetime.date(2019, 10, 14)
    R = 4

    Nx,Ny = 40,30
    xBins = np.logspace(0, 6, Nx+1)
    yBins = np.linspace(-200,50,Ny+1)
    # fBins = np.array([11, 12, 14, 15, 16, 17, 18, 19, 20])
    fBins = np.array(range(5,21))

    fluxMean_log, fluxStd_log = fit_fluxes(date_start, date_end, R, xBins, yBins, fBins)
    return fluxMean_log


def mainfun():
    date_start = datetime.date(2013, 9, 7)
    date_end = datetime.date(2019, 10, 14)
    R = 3

    Nx,Ny = 40,30
    xBins = np.logspace(0, 6, Nx+1)
    yBins = np.linspace(-200,50,Ny+1)
    # fBins = np.array([11, 12, 14, 15, 16, 17, 18, 19, 20])
    fBins = np.array(range(5,21))

    fluxMean_log, fluxStd_log = fit_fluxes(date_start, date_end, R, xBins, yBins, fBins)
    energy = get_energy(isNotNaN=False)
    energy = energy[fBins]

    # fig, ax = plt.subplots(1,1)
    # X,Y = np.meshgrid(energy,yBins)
    # ph = ax.pcolor(X,Y,fluxMean_log.T,shading='flat')
    # fig.colorbar(ph, ax=ax)
    # ax.set_xscale('log')
    # plt.show()


    fig, ax = plt.subplots(1,1)
    X,Y = np.meshgrid(energy,yBins)
    ph = ax.pcolor(X,Y,fluxMean_log.T,shading='flat')
    fig.colorbar(ph, ax=ax)
    ax.set_xscale('log')
    plt.show()



    # popt  = np.zeros_like

    # for fBin in fBins:
    #     binsN = load_binsN(date_start, date_end, R, xBins, yBins, fBin, 'SymH', 'fluxes1')
    #     popt = fit_fluxes_onesE(xBins, binsN)

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(yBins[1:], popt[:,0])
    # plt.show()


# mainfun()
# a = 1