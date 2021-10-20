import pickle
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyspedas


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

def add_bins(xStr, yStr, zStr, binsN, binsData, out_data, xBins, yBins, zBins):
    Nx,Ny = len(xBins)-1,len(yBins)-1
    x0 = out_data.get(xStr)
    y0 = out_data.get(yStr)

    for k in range(len(zBins)):
        z = out_data.get(zStr)[:,zBins[k]]
        inds0 = np.where(z>0)
        x = x0[inds0]
        y = y0[inds0]
        z = z[inds0]        
        for i in range(Nx):
            for j in range(Ny):
                indsX = np.logical_and(xBins[i] < x, x <= xBins[i+1])
                indsY = np.logical_and(yBins[j] < y, y <= yBins[j+1])
                inds  = np.logical_and(indsX, indsY)
                binsN[k,j,i]    += np.sum(inds)
                binsData[k,j,i] += np.sum(z[inds])
    return binsN, binsData

def load_bins_data(date_start, date_end, R, xBins, yBins, zBins, xNames, yNames, zNames, N_load_step = 60, N_load_step_part = 20):
    ## all shapes of Bins must be equal
    if(len(xNames) != len(yNames) or len(xNames) != len(zNames)):
        raise Exception('List of names lenght problem')
    else:
        Nnames = len(xNames)

    Nx,Ny,Nz = len(xBins[0])-1,len(yBins[0])-1,len(zBins)
    binsN    = np.zeros((Nnames, Nz,Ny,Nx),dtype='int32')
    binsData = np.zeros((Nnames, Nz,Ny,Nx),dtype='float64')

    date_current = date_start
    date_step = datetime.timedelta(days=N_load_step)

    while(date_current < date_end):
        date_start_p = date_current
        date_current += date_step
        date_end_p = min(date_current, date_end)
        out_data = load_outdata_partial(date_start_p, date_end_p, R, days_step=N_load_step_part)
        for i in range(Nnames):
            binsN[i,:,:,:], binsData[i,:,:,:] = add_bins(xNames[i], yNames[i], zNames[i], binsN[i,:,:,:], binsData[i,:,:,:], out_data, xBins[i], yBins[i], zBins)

    zData = np.zeros_like(binsData)
    inds = np.where(binsN>0)
    zData[inds] = binsData[inds] / binsN[inds]
    return zData

def plot2dHyst(fig,ax, xBins, yBins, zData):
    X,Y = np.meshgrid(xBins,yBins)
    ph = ax.pcolor(X,Y,zData,shading='flat',norm=matplotlib.colors.LogNorm())
    fig.colorbar(ph, ax=ax)

def makeplot_fluxes(date_start, date_end, R):
    Nx,Ny = 40,30
    # xBins = [np.linspace(45,70,Nx+1), np.linspace(-70,-45,Nx+1)]
    xBins = [np.linspace(1,10,Nx+1), np.linspace(1,10,Nx+1)]
    yBins = [np.linspace(-200,50,Ny+1), np.linspace(-200,50,Ny+1)]
    zBins = [3,9,10,13,14]
    xNames = ['Lshell', 'Lshell']
    yNames = ['SymH', 'SymH']
    zNames = ['fluxes1', 'fluxes2']
    N_load_step, N_load_step_part = 15, 5
    zData = load_bins_data(date_start, date_end, R, xBins, yBins, zBins, xNames, yNames, zNames, N_load_step, N_load_step_part)
    
    zData01_1 = zData[0,0,:,:]
    zData05_1 = (zData[0,1,:,:] + zData[0,2,:,:])/2
    zData10_1 = (zData[0,3,:,:] + zData[0,4,:,:])/2

    zData01_2 = zData[1,0,:,:]
    zData05_2 = (zData[1,1,:,:] + zData[1,2,:,:])/2
    zData10_2 = (zData[1,3,:,:] + zData[1,4,:,:])/2

    fig, axs = plt.subplots(3,2)
    plot2dHyst(fig,axs[0,0], xBins[0], yBins[0], zData01_1)
    plot2dHyst(fig,axs[1,0], xBins[0], yBins[0], zData05_1)
    plot2dHyst(fig,axs[2,0], xBins[0], yBins[0], zData10_1)
    plot2dHyst(fig,axs[0,1], xBins[1], yBins[1], zData01_2)
    plot2dHyst(fig,axs[1,1], xBins[1], yBins[1], zData05_2)
    plot2dHyst(fig,axs[2,1], xBins[1], yBins[1], zData10_2)
    
    axs[0,0].set_title('North fluxes')
    axs[0,0].set_ylabel('E = 100keV\n%s' % yNames[0])
    axs[1,0].set_ylabel('E = 500keV\n%s' % yNames[0])
    axs[2,0].set_ylabel('E = 1MeV\n%s' % yNames[0])
    axs[2,0].set_xlabel(xNames[0])
    axs[0,0].set_xticklabels([])
    axs[1,0].set_xticklabels([])

    axs[0,1].set_title('South fluxes')
    axs[2,1].set_xlabel(xNames[1])
    axs[0,1].set_xticklabels([])
    axs[0,1].set_yticklabels([])
    axs[1,1].set_xticklabels([])
    axs[1,1].set_yticklabels([])
    axs[2,1].set_yticklabels([])

    plt.subplots_adjust(left=0.09, bottom=0.05, right=0.98, top=0.95, wspace=0, hspace=0.1)

    plt.show()


date_start = datetime.date(2015, 1, 1)
date_end = datetime.date(2017, 1, 1)
R = 3
makeplot_fluxes(date_start, date_end, R)
