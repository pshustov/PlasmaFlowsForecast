import pickle
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyspedas
import os
import hashlib


def dict_conc(A, B):
    if(A):
        attr_names_A = list(A.keys()) 
        attr_names_B = list(B.keys()) 
        if not attr_names_A == attr_names_B:
            raise Exception('Different set of keys in dictionaries')
        attr_names_nonNParray = list()
        for attr_name in attr_names_A:
            if isinstance(A.get(attr_name),np.ndarray):
                exec("A.update({'%s': np.concatenate((A.get('%s'), B.get('%s')), axis=0)})" % (attr_name, attr_name, attr_name))
            else:
                attr_names_nonNParray.append(attr_name)
        if(attr_names_nonNParray != list()):
            print(attr_names_nonNParray)
            raise Exception('Unknown class attribute')
    else:
        A = B
    return A


def load_outdata_partial(date_start, date_end, R, days_step=30):
    dirname = './outdata'
    date_current = date_start
    one_month = datetime.timedelta(days=days_step)
    out_data = {}
    while(date_current < date_end):
        date_start_p = date_current
        date_current += one_month
        date_end_p = min(date_current, date_end)
        out_data = dict_conc(out_data,load_outdata(date_start_p, date_end_p, R))
    return out_data


def load_outdata(date_start, date_end, R):
    dirname = './outdata'
    date_current = date_start
    one_day = datetime.timedelta(days=1)
    out_data = {}
    while(date_current < date_end):
        filename = dirname + '/outdata_' + str(R) + '_' + date_current.strftime('%Y-%m-%d.dat')
        with open(filename, 'rb') as fileData:
            out_data = dict_conc(out_data,pickle.load(fileData))
        date_current += one_day
    return out_data


def add_bins(yName, zName, binsN, out_data, xBins, yBins, fBins):
    Nx = len(xBins)-1
    Ny = len(yBins)-1
    y0 = out_data.get(yName)

    for fBin in fBins:
        f = out_data.get(zName)[:,fBin]
        inds0 = np.where(f>0)
        y = y0[inds0]
        f = f[inds0]     
        for i in range(Nx):   
            for j in range(Ny):
                indsX = np.logical_and(xBins[i] < f, f <= xBins[i+1])
                indsY = np.logical_and(yBins[j] < y, y <= yBins[j+1])
                inds  = np.logical_and(indsX, indsY)
                binsN[j,i] += np.sum(inds)
    return binsN


def load_bins_data(date_start, date_end, R, xBins, yBins, fBin, yName, zName, N_load_step = 60, N_load_step_part = 20):
    Nx,Ny    = len(xBins)-1,len(yBins)-1
    
    dirname = './data_distrib'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    hashMD5 = hashlib.md5()
    
    hashMD5.update(date_start.strftime('%Y%m%d').encode('utf-8'))
    hashMD5.update(date_end.strftime('%Y%m%d').encode('utf-8'))
    hashMD5.update(str(R).encode('utf-8'))
    hashMD5.update(xBins.astype('float32').tobytes())
    hashMD5.update(yBins.astype('float32').tobytes())
    hashMD5.update(np.ndarray(fBin).astype('int32').tobytes())
    hashMD5.update(yName.encode('utf-8'))
    hashMD5.update(zName.encode('utf-8'))


    filename = dirname+'/data_distrib_'+hashMD5.hexdigest()+'.dat'
    
    if not os.path.isfile(filename):
        binsN    = np.zeros((Ny,Nx),dtype='int32')
        date_current = date_start
        date_step = datetime.timedelta(days=N_load_step)

        while(date_current < date_end):
            date_start_p = date_current
            date_current += date_step
            date_end_p = min(date_current, date_end)
            out_data = load_outdata_partial(date_start_p, date_end_p, R, days_step=N_load_step_part)        
            binsN = add_bins(yName, zName, binsN, out_data, xBins, yBins, fBin)
            print(date_current)
        
        with open(filename, 'wb') as fileData:
            pickle.dump(binsN,fileData,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileData:
            binsN = pickle.load(fileData)

    binsData = (binsN.T / np.sum(binsN,axis=1)).T
    return binsData


def plot2dHyst(fig,ax, xBins, yBins, zData):
    X,Y = np.meshgrid(xBins,yBins)
    ph = ax.pcolor(X,Y,zData,shading='flat',norm=matplotlib.colors.LogNorm())
    fig.colorbar(ph, ax=ax)


def makeplot_fluxes_distrib_onDST(fig, ax, date_start, date_end, R, yName, zName, fBin, Nx=40, Ny=40):
    xBins = np.logspace(0, 6, Nx+1)
    yBins = np.linspace(-200,50,Ny+1)

    N_load_step, N_load_step_part = 60, 20
    zData = load_bins_data(date_start, date_end, R, xBins, yBins, fBin, yName, zName, N_load_step, N_load_step_part)

    plot2dHyst(fig,ax, xBins, yBins, zData)
    ax.set_xscale('log')


date_start = datetime.date(2013, 1, 1)
date_end = datetime.date(2017, 10, 14)
R = 4

Nx,Ny = 40,40

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