import numpy as np
import pyspedas
import PyGeopack as gp
import datetime
import pickle
import os
import hashlib
import kpindex
import pyomnidata


from matplotlib import pyplot as plt
from Globals import *
from PyGeopack.TraceField import TraceField
from scipy.interpolate import interp1d
from math import ceil
from scipy.optimize.minpack import curve_fit
from scipy.special import gamma

# gp.UpdateParameters(SkipWParameters=True)

class getOmni():
    def __init__(self, date=None):
        if date is None:
            self.year = None
            self.dataOmni = None
            self.omni_time = None
        else:
            self.year = date.year
            self.dataOmni = pyomnidata.ReadOMNI(date.year)
            utc = self.dataOmni['utc']
            self.omni_time = np.zeros_like(utc)
            for i in range(len(utc)):
                self.omni_time[i] = (datetime.datetime(year=1950,month=1,day=1) + datetime.timedelta(hours=utc[i])).timestamp()
    
    def update(self, date):
        if(self.year != date.year):
            self.__init__(date)

    def get_data(self, date):
        self.update(date)
        return self.omni_time, self.dataOmni
            

class getKp():
    def __init__(self):
        kpInd = kpindex.GetKp()
        dates = kpInd['Date']
        ut0 = kpInd['ut0'].astype(np.float64)
        self.kp_time = np.zeros_like(dates)
        for i in range(len(dates)):
            self.kp_time[i] = (datetime.datetime.strptime(str(dates[i]),'%Y%m%d') + datetime.timedelta(hours=ut0[i])).timestamp()
        self.kp_data = kpInd['Kp']
    
    def get_data(self):
        return self.kp_time, self.kp_data


Kp = getKp()
omni = getOmni()


class TraceFieldMy(TraceField):
    def __init__(self, times = None, positions = None, Model='T96', CoordIn = 'GSM', CoordOut = 'GSM', 
				alt=100.0, MaxLen=1000, DSMax=1.0,FlattenSingleTraces=True,Verbose=False,OutDtype='float64',**kwargs):
        if times is not None:
            TraceField.__init__(self,*positions.T, *times_to_dateUt(times), Model, CoordIn, CoordOut, alt, MaxLen, DSMax,FlattenSingleTraces,Verbose,OutDtype,**kwargs)
            self.time = times
            self.pos  = positions
            self.nstep = self.nstep.astype('int32')

    def append(self, other):
        if(self.__dict__):
            attr_names = list(self.__dict__.keys()) 
            attr_names_nonNParray = list()
            for attr_name in attr_names:
                if isinstance(self.__dict__.get(attr_name),np.ndarray):
                    exec('self.%s = np.concatenate((self.%s, other.%s), axis=0)' % (attr_name, attr_name, attr_name))
                else:
                    attr_names_nonNParray.append(attr_name)
            if(attr_names_nonNParray != list('n')):
                raise Exception('Unknown class attribute')
            else:
                self.n += other.n
        else:                
            self.__dict__ = other.__dict__


def times_to_dateUt(times):
        pos_time_DT = [ datetime.datetime.utcfromtimestamp(times[i]) for i in range(len(times))]
        Date        = np.array([int(pos_time_DT[i].strftime('%Y%m%d')) for i in range(len(times))])
        ut          = np.array([(pos_time_DT[i] - datetime.datetime(pos_time_DT[i].year,pos_time_DT[i].month,pos_time_DT[i].day)).total_seconds() for i in range(len(times))]) / 3600
        return Date,ut


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


def get_time_inds(pos_time, trange = None, t_step = 1, R = None):
    n_step = int(ceil(t_step.total_seconds()/np.mean(np.diff(pos_time))))
    inds = np.full(np.size(pos_time), False, dtype=bool)
    if R is None:
        inds[::n_step] = True
    else:
        R = R / 3
        if trange is not None:
            inds = pos_time > pyspedas.time_float(trange[0])
            i = np.argwhere(pos_time > pyspedas.time_float(trange[0]))[0]
        else:
            i = 0
        while(i < len(inds)):
            inds[i] = True
            i += ceil(n_step * R[i])    
    inds[-1] = True

    if trange is None:
        return inds
    else:
        return (pos_time >= pyspedas.time_float(trange[0])) & (pos_time < pyspedas.time_float(trange[1])) & inds


def interp1(x,y,xp):
    f = interp1d(x,y,bounds_error=False,fill_value=np.nan)
    return f(xp)


def get_TraceMy(date, pos_time=None, pos_position=None):
    t_step   = datetime.timedelta(minutes=2)

    dirname = database+'rbsp_magtraces/'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    filename        = dirname + date.strftime('magtrace%Y-%m-%d.dat')
    if not os.path.isfile(filename):
        if (pos_time is None) and (pos_position is None):
            dateCurrentStr  = date.strftime('%Y-%m-%d/%H:%M:%S')
            dateNextStr     = (date+datetime.timedelta(days=1)).strftime('%Y-%m-%d/%H:%M:%S')
            emfisis_vars    = pyspedas.rbsp.emfisis([dateCurrentStr, dateNextStr], probe='a', level='l3', notplot=True, coord='gsm')

            pos_time     = emfisis_vars.get('coordinates').get('x')
            pos_position = emfisis_vars.get('coordinates').get('y') / 6371.2 

        wherenan = np.where(np.logical_not(np.isnan(pos_position[:,1])))
        pos_time = pos_time[wherenan]
        pos_position = pos_position[wherenan,:][0]

        pos_R = np.sqrt(pos_position[:,0]**2 + pos_position[:,1]**2 + pos_position[:,2]**2)
        inds         = get_time_inds(pos_time, t_step=t_step, R=pos_R)
        pos_time     = pos_time[inds]
        pos_position = pos_position[inds,:]

        Trace = TraceFieldMy(pos_time, pos_position,Verbose=True)
        with open(filename, 'wb') as fileTrace:
            pickle.dump(Trace, fileTrace,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileTrace:
            Trace = pickle.load(fileTrace)
    return Trace


def get_coord_intersect(out_time, R, Trace):    
    trace_time      = Trace.time
    trace_x         = Trace.x
    trace_y         = Trace.y
    trace_z         = Trace.z
    trace_R         = np.sqrt(trace_x**2 + trace_y**2 + trace_z**2)
    inds            = np.nanargmax(trace_R, axis=1)
    n = len(trace_time)
    
    x1, x2 = np.zeros(n), np.zeros(n)
    y1, y2 = np.zeros(n), np.zeros(n)
    z1, z2 = np.zeros(n), np.zeros(n)
    for i in range(n):
        x1[i] = interp1(trace_R[i,:inds[i]+1],trace_x[i,:inds[i]+1], R)
        x2[i] = interp1(trace_R[i,inds[i]:],trace_x[i,inds[i]:], R)
        y1[i] = interp1(trace_R[i,:inds[i]+1],trace_y[i,:inds[i]+1], R)
        y2[i] = interp1(trace_R[i,inds[i]:],trace_y[i,inds[i]:], R)
        z1[i] = interp1(trace_R[i,:inds[i]+1],trace_z[i,:inds[i]+1], R)
        z2[i] = interp1(trace_R[i,inds[i]:],trace_z[i,inds[i]:], R)

    X1 = interp1(trace_time,x1,out_time)
    X2 = interp1(trace_time,x2,out_time)
    Y1 = interp1(trace_time,y1,out_time)
    Y2 = interp1(trace_time,y2,out_time)
    Z1 = interp1(trace_time,z1,out_time)
    Z2 = interp1(trace_time,z2,out_time)
    return X1,Y1,Z1, X2,Y2,Z2


def get_fluxes_integrated(mageis_vars, X,Y,Z):
    FEDU_t = mageis_vars.get('FEDU').get('x')
    fluxes = mageis_vars.get('FEDU').get('y')
    alpha0 = mageis_vars.get('FEDU').get('v1') 
    energy = mageis_vars.get('FEDU').get('v2') 
    position = mageis_vars.get('Position').get('y')
    if(np.nanmean(np.abs(position)) > 40):
        position = position / 6371.2
    
    Bx0,By0,B0z = gp.ModelField(*position.T,*times_to_dateUt(FEDU_t),Model='T96',CoordIn='GSM',CoordOut='GSM')
    Bx1,By1,Bz1 = gp.ModelField(X,Y,Z,*times_to_dateUt(FEDU_t),Model='T96',CoordIn='GSM',CoordOut='GSM')
    B0 = np.sqrt(Bx0**2 + By0**2 + B0z**2)
    B1 = np.sqrt(Bx1**2 + By1**2 + Bz1**2)

    alpha = np.arcsin(np.outer(np.sqrt(B1 / B0), np.sin(alpha0 / 180 * np.pi))) * 180 / np.pi
    alpha[:,np.where(alpha0>90)] = 180 - alpha[:,np.where(alpha0>90)]

    fluxes[np.where(np.isnan(alpha))] = np.nan
    fluxes = np.nansum(fluxes,axis=1)
    return fluxes #, alpha, energy


def interp1_param(x, y, xp):
    inds =np.logical_not(np.isnan(y))
    x = x[inds]
    y = y[inds]
    return interp1(x,y,xp)


def get_params(out_time, date):
    omni_time, dataOmni = omni.get_data(date)
    kp_time, kp_data    = Kp.get_data()

    Pdyn = interp1_param(omni_time, dataOmni['FlowPressure'], out_time)
    SymH = interp1_param(omni_time, dataOmni['SymH'],         out_time)
    By   = interp1_param(omni_time, dataOmni['ByGSM'],        out_time)
    Bz   = interp1_param(omni_time, dataOmni['BzGSM'],        out_time)  
    kp   = interp1_param(kp_time, kp_data, out_time)
    return Pdyn, SymH, By, Bz, kp


def get_projectionAt_oneDay(date, R):
    """ 
    Get data projected at R for one day

    R: double -- radius of a sphere in RE at which projection should be

    date: datetime.date  --  date for projection 
    """
    dirname = database+'outdata/'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    filename = dirname + 'outdata_' + str(R) + '_' + date.strftime('%Y-%m-%d.dat')
    
    if not os.path.isfile(filename):
        date_current_str  = date.strftime('%Y-%m-%d/%H:%M:%S')
        date_next_str     = (date+datetime.timedelta(days=1)).strftime('%Y-%m-%d/%H:%M:%S')  
        mageis_vars = pyspedas.rbsp.mageis([date_current_str, date_next_str], probe='a', level='l3', rel='rel04', notplot=True)
        out_time = mageis_vars.get('FEDU').get('x')
        pos_time = mageis_vars.get('Position').get('x')
        pos_position = mageis_vars.get('Position').get('y')
        if(np.nanmean(np.abs(pos_position)) > 40):
            pos_position = pos_position / 6371.2
        

        Trace   = get_TraceMy(date,pos_time,pos_position)
        X1,Y1,Z1, X2,Y2,Z2 = get_coord_intersect(out_time, R, Trace)
        
        out_data = {}
        
        Pdyn, SymH, By, Bz, kp = get_params(out_time, date)
        out_data.update({'MLAT1':   interp1(Trace.time, Trace.MlatN, out_time)})
        out_data.update({'MLAT2':   interp1(Trace.time, Trace.MlatS, out_time)}) 
        out_data.update({'Lshell':  interp1(Trace.time, Trace.Lshell, out_time)})
        out_data.update({'fluxes1': get_fluxes_integrated(mageis_vars, X1,Y1,Z1)})
        out_data.update({'fluxes2': get_fluxes_integrated(mageis_vars, X2,Y2,Z2)})  
        out_data.update({'Pdyn':    Pdyn})
        out_data.update({'SymH':    SymH})
        out_data.update({'By':      By})
        out_data.update({'Bz':      Bz})
        out_data.update({'kp':      kp})
        
        with open(filename, 'wb') as fileData:
            pickle.dump(out_data,fileData,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileData:
            out_data = pickle.load(fileData)
    return out_data


def get_projectionAt(R, date_start=datetime.date(2012, 9, 7), date_end=datetime.date(2019, 10, 14),  days_step=30):
    """ 
    Get data projected at R from date_start to date_end

    R: double -- radius of a sphere in RE at which projection should be

    date_start: datetime.date  --  start date
    
    date_end: datetime.date  --  end date 

    days_step: uint  --  summing of dictionary routine interval in days
    """
    date_current = date_start
    one_day = datetime.timedelta(days=1)
    date_step = datetime.timedelta(days=days_step)
    out_data = {}
    
    while(date_current < date_end):
        date_end_p = min(date_current+date_step, date_end)
        out_data_p = {}
        while(date_current < date_end_p):
            out_data_p = dict_conc(out_data_p, get_projectionAt_oneDay(date_current, R))
            date_current += one_day
        out_data = dict_conc(out_data,out_data_p)

    return out_data


def add_bins(bins, data, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N):
    """
    Function calculate distribution number of bins, based on "paramName" and fluxes and just add tham fo bins

    bins: ndarray  --  array which contain bin numbers, may be nonzero

    data: ndarray  --  dictionary of input data 

    paramName: str  --  name of parameter

    paramBins_b: ndaray  --  vector of boundaries for parameter bins 

    fluxesName: str  --  name of fluxes ("fluxes1" for north and "fluxes2" for south)

    fluxesBins_b: ndaray  --  vector of boundaries for fluxes bins

    energyBins_N: iterable  --  array of energy bins with int numbers
    """
    Nx = len(fluxesBins_b)-1
    Ny = len(paramBins_b)-1
    y0 = data.get(paramName)

    for eBin in energyBins_N:
        f = data.get(fluxesName)[:,eBin]
        inds0 = np.where(f>0)
        y = y0[inds0]
        f = f[inds0]     
        for i in range(Nx):   
            for j in range(Ny):
                indsX = np.logical_and(fluxesBins_b[i] < f, f <= fluxesBins_b[i+1])
                indsY = np.logical_and(paramBins_b[j] < y, y <= paramBins_b[j+1])
                inds  = np.logical_and(indsX, indsY)
                bins[j,i] += np.sum(inds)
    return bins


def load_binsN(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N, N_load_step = 60, N_load_step_part = 20):
    """
    Function return bins array, with distribution on fluxes and "paramName"

    R: double  --  radius of a projection sphere 

    date_start: datetime.date  --  start date of a projection data

    date_end: datetime.date  --  end date of a projection data

    paramName: str  --  name of parameter

    paramBins_b: ndaray  --  vector of boundaries for parameter bins 

    fluxesName: str  --  name of fluxes ("fluxes1" for north and "fluxes2" for south)

    fluxesBins_b: ndaray  --  vector of boundaries for fluxes bins

    energyBins_N: iterable  --  array of energy bins with int numbers
    """
    Nx,Ny    = len(fluxesBins_b)-1,len(paramBins_b)-1
    
    dirname = database+'data_binsdistrib/'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    hashMD5 = hashlib.md5()
    
    hashMD5.update(np.array(R).astype('float16').tobytes())
    hashMD5.update(date_start.strftime('%Y%m%d').encode('utf-8'))
    hashMD5.update(date_end.strftime('%Y%m%d').encode('utf-8'))
    hashMD5.update(paramName.encode('utf-8'))
    hashMD5.update(paramBins_b.astype('float16').tobytes())
    hashMD5.update(fluxesName.encode('utf-8'))
    hashMD5.update(fluxesBins_b.astype('float16').tobytes())
    hashMD5.update(np.array(energyBins_N).astype('int16').tobytes())

    # energyBins_all = np.array([i for i in range(len(get_energy_rbsp_mageis()))])
    filename = dirname+'binsdistrib_'+hashMD5.hexdigest()+'.dat'
    
    if not os.path.isfile(filename):
        binsN    = np.zeros((Ny,Nx),dtype='int32')
        date_current = date_start
        date_step = datetime.timedelta(days=N_load_step)

        while(date_current < date_end):
            date_start_p = date_current
            date_current += date_step
            date_end_p = min(date_current, date_end)
            out_data = get_projectionAt(R, date_start_p, date_end_p, days_step=N_load_step_part)        
            binsN = add_bins(binsN, out_data, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N)
            print("Bins calculated:" + date_current.strftime('%Y-%m-%d'))
        
        with open(filename, 'wb') as fileData:
            pickle.dump(binsN,fileData,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileData:
            binsN = pickle.load(fileData)

    return binsN


# def logfun(x, x0, k, a0):
#     (a0*k+1) * np.log(x/x0) + (-k) * np.log(1 + x/x0) - np.log(gamma(a0*k)*gamma(k - a0*k))/gamma(k))


def fit_fluxes_gauss(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N, is_save_figures=False, is_save_cvs = False):
    """
    Function fits bins by normal gauss, return fluxMean_log and fluxStd_log

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


    fluxMean_log = np.zeros((len(energyBins_N),len(paramBins_b)-1))
    fluxStd_log  = np.zeros((len(energyBins_N),len(paramBins_b)-1))
    
    xFluxes = fluxesBins_b[0:-1] + np.diff(fluxesBins_b) / 2

    fun  = lambda x, x0, k, a0: (x/x0)**(a0*k+1) * (1 + (x/x0))**(-k) / (gamma(a0*k)*gamma(k - a0*k)/gamma(k))
    logfun = lambda x, x0, k, a0: np.log(fun(x, x0, k, a0))

    poptM = np.zeros((len(energyBins_N),len(paramBins_b)-1, 3))
    energy = get_energy_rbsp_mageis()

    if is_save_figures:
        fig = plt.figure(figsize=[6.5,19.5])
        fig_dirname = database+'figures\\step_fit\\'
        if not os.path.isdir(fig_dirname):
            os.makedirs(fig_dirname)
        paramBinsNames = [('%.0f..%.0f'%(par1,par2)) for par1,par2 in zip(paramBins_b[:-1],paramBins_b[1:])]
        paramBinsNames.append(paramName)


    if is_save_cvs:
        cvsdirname = '.\\cvs\\'
        if not os.path.isdir(cvsdirname):
            os.makedirs(cvsdirname)
        np.savetxt(cvsdirname+'energy.cvs', energy[energyBins_N], fmt='%.8e', delimiter=',')
        np.savetxt(cvsdirname+'fluxesBins_b.cvs', fluxesBins_b, fmt='%.8e', delimiter=',')


    for i in range(len(energyBins_N)):
        binsN    = load_binsN(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, np.array([energyBins_N[i]]))
        
        if is_save_cvs:
            np.savetxt(cvsdirname+('binsN%d.cvs'%(i)), binsN, fmt='%.8e', delimiter=',')

        norm     = np.sum(binsN,axis=1)
        norm[np.where(norm<5)] = -1
        binsData = (binsN.T / norm).T
        
        for j in range(binsData.shape[0]):
            if(norm[j]>0):
                indsXb = np.argwhere(binsN[j,:]>2)
                p0 = [ xFluxes[np.argmax(binsData[j,:])], 10, 0.5]
                bounds = ((1e-6, 0, 0), (1e6, 100, 1))
                yData = np.log(binsData[j,:])
                indsY = np.where(np.logical_and(np.logical_and(~np.isinf(yData), xFluxes>0.05), binsData[j,:]>1000/norm[j]))
                sigma = 1 / np.sqrt(binsN[j,:])
                popt, pcov = curve_fit(logfun, xFluxes[indsY], yData[indsY], p0=p0, bounds=bounds)
            else:
                popt = (0,0)
            poptM[i,j,:] = popt
        
        if is_save_figures:
            fig_name = ('figure_distrib_fit_E=%gkeV.png' % energy[energyBins_N[i]])
            ax = fig.add_axes([0.13, 0.02, 0.85, 0.94])
            for j in range(binsData.shape[0]):
                norm_binsdata = np.max(binsData[j,:]) / 0.95
                norm2 = 3
                stepY = np.log10(binsData[j,:]/norm_binsdata) / norm2 + j + 1
                dataY = np.log10(fun(xFluxes, *poptM[i,j,:])/norm_binsdata)
                plotY = dataY / norm2 + j + 1
                indsY = np.where(dataY>-1*norm2-1)

                ax.step(fluxesBins_b[1:], stepY, where='pre', color='tab:blue')
                ax.plot(xFluxes[indsY], plotY[indsY], color='tab:orange')
                ax.axhline(j, color='gray', linewidth=0.5)
                ax.text(fluxesBins_b[-1], j, ('%d'%norm[j]), ha='right', va='bottom')
                ax.set_yticks(np.arange(len(paramBinsNames)))
                ax.set_yticklabels(paramBinsNames)
                ax.set_xscale('log')
            ax.set_title('%g keV' % energy[energyBins_N[i]])
            ax.set_ylim([0,binsData.shape[0]+0.5])
            ax.set_xlim([np.min(fluxesBins_b), np.max(fluxesBins_b)])
            fig.savefig(fig_dirname + fig_name)
            fig.clf()

    if is_save_figures:
        plt.close(fig)

    return fluxMean_log, fluxStd_log, paramBins_b


def get_equinumeros_subests_on_param(paramName, numOfBins, R, date_start=datetime.date(2012, 9, 7), date_end=datetime.date(2019, 10, 14)):
    """
    Return equinumeros distributed bins on nonempty data from date_start to date_end on R sphere based on paramName. 

    paramName: str  --  name of param from projected data on which bins should be constracted 

    numOfBins: uint  --  number of bins to return

    R: double  --  radius of a projection sphere 

    date_start: datetime.date  --  start date of a projection data

    date_end: datetime.date  --  end date of a projection data
    """
    date_current = date_start
    date_step = datetime.timedelta(days=60)

    dirname = database+'equibins/'
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
            out_data = get_projectionAt(R, date_start_p, date_end_p, days_step=20)
            param_data = out_data.get(paramName)
            fluxes_data = out_data.get('fluxes1')
            val = np.concatenate((val,param_data[np.any(fluxes_data > 0,axis=1)]))
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


def get_energy_rbsp_mageis():
    """
    Return energy bins from rbsp mageis.
    """

    dirname = database
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    
    filename = database + 'rbsp_mageis_energies.dat'
    
    if not os.path.isfile(filename):
        mageis_vars = pyspedas.rbsp.mageis(['2015-01-01/00:00:00', '2015-01-01/01:00:00'], probe='a', level='l3', rel='rel04', notplot=True)
        energy = mageis_vars.get('FEDU').get('v2') 
        with open(filename, 'wb') as fileData:
            pickle.dump(energy,fileData,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileData:
            energy = pickle.load(fileData)
    return energy
