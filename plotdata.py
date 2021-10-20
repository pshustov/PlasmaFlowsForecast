import pickle
import datetime
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as tri
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

def load_outdata_partial(date_start, date_end, R):
    dirname = './outdata'
    date_current = date_start
    one_month = datetime.timedelta(days=30)
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

def plot2dHyst(fig,ax,x,y,z,Nx,Ny):
    inds = np.where(z>0)
    print(np.argwhere(z>0).shape[0]/x.shape[0])
    x = x[inds]
    y = y[inds]
    z = z[inds]
    xBins = np.linspace(np.min(x),np.max(x),Nx+1)
    yBins = np.linspace(np.min(y),np.max(y),Ny+1)
    bins = np.zeros((Ny,Nx))
    for i in range(Nx):
        for j in range(Ny):
            indsX = np.logical_and(xBins[i] < x, x <= xBins[i+1])
            indsY = np.logical_and(yBins[j] < y, y <= yBins[j+1])
            bins[j,i] = np.mean(z[np.where(np.logical_and(indsX,indsY))])
    X,Y = np.meshgrid(xBins,yBins)
    ph = ax.pcolor(X,Y,bins,shading='flat',norm=matplotlib.colors.LogNorm())
    fig.colorbar(ph, ax=ax)

date_start = datetime.date(2013, 1, 1)
date_end = datetime.date(2013, 1, 20)
R = 3
out_data = load_outdata_partial(date_start, date_end, R)


# mageis_vars = pyspedas.rbsp.mageis([date_start.strftime('%Y-%m-%d/%H:%M:%S'), (date_start+datetime.timedelta(days=1)).strftime('%Y-%m-%d/%H:%M:%S')  ], probe='a', level='l3', rel='rel04', notplot=True)
# energy = mageis_vars.get('FEDU').get('v2') 
# print(energy[ii])

MLAT1 = out_data.get('MLAT1')
MLAT2 = out_data.get('MLAT2')
SymH  = out_data.get('SymH')
z01_1 = out_data.get('fluxes1')[:,3]
z05_1 = (out_data.get('fluxes1')[:,9] + out_data.get('fluxes1')[:,10]) / 2
z10_1 = (out_data.get('fluxes1')[:,13] + out_data.get('fluxes1')[:,14]) / 2
z01_2 = out_data.get('fluxes2')[:,3]
z05_2 = (out_data.get('fluxes2')[:,9] + out_data.get('fluxes2')[:,10]) / 2
z10_2 = (out_data.get('fluxes2')[:,13] + out_data.get('fluxes2')[:,14]) / 2

Nx,Ny = 60,80
fig, axs = plt.subplots(3,2)
plot2dHyst(fig,axs[0,0],MLAT1,SymH,z01_1,Nx,Ny)
plot2dHyst(fig,axs[1,0],MLAT1,SymH,z05_1,Nx,Ny)
plot2dHyst(fig,axs[2,0],MLAT1,SymH,z10_1,Nx,Ny)
plot2dHyst(fig,axs[0,1],MLAT1,SymH,z01_2,Nx,Ny)
plot2dHyst(fig,axs[1,1],MLAT1,SymH,z05_2,Nx,Ny)
plot2dHyst(fig,axs[2,1],MLAT1,SymH,z10_2,Nx,Ny)
plt.show()
