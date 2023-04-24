import numpy as np
import datetime
import pickle
import os
import hashlib
import dataFlux


from matplotlib import pyplot as plt
from Globals import *
from scipy.optimize.minpack import curve_fit
from scipy.special import gamma



# def logfun(x, x0, k, a0):
#     (a0*k+1) * np.log(x/x0) + (-k) * np.log(1 + x/x0) - np.log(gamma(a0*k)*gamma(k - a0*k))/gamma(k))


def fit_fluxes(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N, is_save_figures=False, is_save_cvs = False):
    """
    Function fits bins

    R: double  --  radius of a projection sphere 

    date_start: datetime.date  --  start date of a projection data

    date_end: datetime.date  --  end date of a projection data

    paramName: str  --  name of parameter

    paramBins_b: ndaray  --  vector of boundaries for parameter bins 

    fluxesName: str  --  name of fluxes ("fluxes1" for north and "fluxes2" for south)

    fluxesBins_b: ndaray  --  vector of boundaries for fluxes bins

    energyBins_N: iterable  --  array of energy bins with int numbers
    """
    if(isinstance(paramBins_b, int)):
        paramBins_b = get_equinumeros_subests_on_param(paramName, paramBins_b, R, date_start, date_end)
    elif(not isinstance(paramBins_b, np.ndarray)):
        raise("Wrong 'paramBins_b'")
    
    fluxes = fluxesBins_b[0:-1] + np.diff(fluxesBins_b) / 2

    numberOfParams = 4
    logfun = lambda x, a, b, p, x0: a - b*x - np.exp(p*(x-x0))

    poptM = np.zeros((len(energyBins_N),len(paramBins_b)-1, numberOfParams))
    energies = dataFlux.get_energy_rbsp_mageis()

    if is_save_figures:
        fig = plt.figure(figsize=[1.5*6.5,1.5*19.5])
        fig_dirname = database+'figures\\step_fit\\'
        if not os.path.isdir(fig_dirname):
            os.makedirs(fig_dirname)
        paramBinsNames = [('%.0f..%.0f'%(par1,par2)) for par1,par2 in zip(paramBins_b[:-1],paramBins_b[1:])]
        paramBinsNames.append(paramName)


    if is_save_cvs:
        cvsdirname = '.\\cvs\\'
        if not os.path.isdir(cvsdirname):
            os.makedirs(cvsdirname)
        np.savetxt(cvsdirname+'energies.cvs', energies[energyBins_N], fmt='%.8e', delimiter=',')
        np.savetxt(cvsdirname+'fluxesBins_b.cvs', fluxesBins_b, fmt='%.8e', delimiter=',')
        np.savetxt(cvsdirname+'paramBins_b.cvs', paramBins_b, fmt='%.8e', delimiter=',')


    for i in range(len(energyBins_N)):
        binsN =1 # = load_binsN(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, np.array([energyBins_N[i]]))
        
        if is_save_cvs:
            np.savetxt(cvsdirname+('binsN%d.cvs'%(i)), binsN, fmt='%.8e', delimiter=',')

        norm = np.sum(binsN,axis=1)
        norm[np.where(norm<1)] = -1
        binsData = (binsN.T / norm).T
        
        for j in range(binsData.shape[0]):
            if(norm[j]>0):
                yData = np.log(binsData[j,:])
                iYmax = np.argmax(binsData[j,:])
                indsData = np.where(np.logical_and(fluxes>fluxes[iYmax], binsN[j,:]>7))
                
                p0 = [ 15, 1, 1, 15]
                bounds = ((-np.inf, 0, 0, -np.inf), (np.inf, np.inf, np.inf, np.inf))
                sigma = 1 / np.sqrt(binsN[j,:])
                
                try:
                    popt, pcov = curve_fit(logfun, np.log(fluxes[indsData]), yData[indsData], p0=p0, bounds=bounds)
                except:
                    popt = (0, 0, 0, 0)
            else:
                popt = (0,0)
            poptM[i,j,:] = popt
        
        if is_save_figures:
            fig_name = ('figure_distrib_fit_E=%gkeV.png' % energies[energyBins_N[i]])
            ax = fig.add_axes([0.13, 0.02, 0.85, 0.94])
            for j in range(binsData.shape[0]):
                norm_binsdata = np.max(binsData[j,:]) / 0.95
                norm2 = 3
                stepY = np.log10(binsData[j,:]/norm_binsdata) / norm2 + j + 1
                dataY = np.log10(np.exp(logfun(np.log(fluxes), *poptM[i,j,:]))/norm_binsdata)
                plotY = dataY / norm2 + j + 1
                indsData = np.where(np.logical_and(dataY>-1*norm2-1, fluxes>=fluxes[np.argmax(stepY)]))

                ax.step(fluxesBins_b[1:], stepY, where='pre', color='tab:blue')
                ax.plot(fluxes[indsData], plotY[indsData], color='tab:orange')
                ax.axhline(j, color='gray', linewidth=0.5)
                ax.text(fluxesBins_b[-1], j, ('%d'%norm[j]), ha='right', va='bottom')
                ax.set_yticks(np.arange(len(paramBinsNames)))
                ax.set_yticklabels(paramBinsNames)
                ax.set_xscale('log')
            ax.set_title('%g keV' % energies[energyBins_N[i]])
            ax.set_ylim([0,binsData.shape[0]+0.5])
            ax.set_xlim([np.min(fluxesBins_b), np.max(fluxesBins_b)])
            fig.savefig(fig_dirname + fig_name)
            fig.clf()

    if is_save_figures:
        plt.close(fig)

    return paramBins_b


def get_equinumeros_subests_on_param(paramName, numOfBins, R, date_start_str="20120907", date_end_str="20191014"):
    """
    Return equinumeros distributed bins on nonempty data from date_start to date_end on R sphere based on paramName. 

    paramName: str  --  name of param from projected data on which bins should be constracted 

    numOfBins: uint  --  number of bins to return

    R: double  --  radius of a projection sphere 

    date_start: datetime.date  --  start date of a projection data

    date_end: datetime.date  --  end date of a projection data
    """

    date_start  = datetime.date(int(date_start_str[0:4]), int(date_start_str[4:6]), int(date_start_str[6:8]))
    date_end    = datetime.date(int(date_end_str[0:4]), int(date_end_str[4:6]), int(date_end_str[6:8]))


    date_current = date_start
    date_step = datetime.timedelta(days=60)

    dirname = database+'equibins\\'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    hashMD5 = hashlib.md5()    
    hashMD5.update(paramName.encode('utf-8'))
    hashMD5.update(np.array(numOfBins).astype('int16').tobytes())
    hashMD5.update(np.array(R).astype('float16').tobytes())
    hashMD5.update(date_start.strftime('%Y%m%d').encode('utf-8'))
    hashMD5.update(date_end.strftime('%Y%m%d').encode('utf-8'))
    filename =  dirname + 'equibins_'+hashMD5.hexdigest()+'.dat'

    if not os.path.isfile(filename):
        val = np.empty(shape=(0,))

        while(date_current < date_end):
            date_start_p = date_current
            date_current += date_step
            date_end_p = min(date_current, date_end)
            out_data =  dataFlux.get_dataflux(R, date_start_p.strftime('%Y%m%d'), date_end_p.strftime('%Y%m%d'), days_step=20)
            val = np.concatenate( (val, out_data.get(paramName)) )
            print('Equinumeros subest: %s' % (date_current))

        val = np.sort(val[~np.isnan(val)])
        N = len(val)
        dN = int(N / numOfBins)
        bins = np.zeros(numOfBins+1)
        for i in range(numOfBins):
            bins[i] = val[dN*i]
        bins[-1] = val[-1]
        with open(filename, 'wb') as fileData:
            pickle.dump(bins,fileData,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileData:
            bins = pickle.load(fileData)
    return bins


get_equinumeros_subests_on_param('SymH', numOfBins=40, R=4, date_start_str="20120907", date_end_str="20130101")