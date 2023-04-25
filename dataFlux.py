import numpy as np
from numpy.lib.index_tricks import AxisConcatenator
import pyspedas
import PyGeopack as gp
import datetime
import pickle
import os
import kpindex
import pyomnidata
import TraceFieldMy as TFm


from matplotlib import pyplot as plt
from Globals import *
from scipy.interpolate import interp1d
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


def dict_conc(A, B):
    if(A):
        if(B):
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
            return A
    else:
        A = B
    return A


def interp1(x,y,xp,defval=np.nan):
    inds =np.logical_not(np.isnan(y))
    x = x[inds]
    y = y[inds]
    f = interp1d(x,y,bounds_error=False,fill_value=defval)
    return f(xp)


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
        y1[i] = interp1(trace_R[i,:inds[i]+1],trace_y[i,:inds[i]+1], R)
        z1[i] = interp1(trace_R[i,:inds[i]+1],trace_z[i,:inds[i]+1], R)

        x2[i] = interp1(trace_R[i,inds[i]:],trace_x[i,inds[i]:], R)
        y2[i] = interp1(trace_R[i,inds[i]:],trace_y[i,inds[i]:], R)
        z2[i] = interp1(trace_R[i,inds[i]:],trace_z[i,inds[i]:], R)

    if(out_time.size):
        pos1 = np.column_stack(( interp1(trace_time,x1,out_time), interp1(trace_time,y1,out_time), interp1(trace_time,z1,out_time) ))
        pos2 = np.column_stack(( interp1(trace_time,x2,out_time), interp1(trace_time,y2,out_time), interp1(trace_time,z2,out_time) ))
    else:
        pos1 = np.array([[], [], []]).T
        pos2 = np.array([[], [], []]).T

    return pos1, pos2


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


def get_coord_GSM(date_start_str, date_end_str, times):
    emfisis_vars = pyspedas.rbsp.emfisis([date_start_str, date_end_str], probe='a', level='l3', notplot=True, coord='gsm')
    emfisis_vars_coord = emfisis_vars.get('coordinates')
    if(emfisis_vars_coord):
        pos_times = np.array(emfisis_vars_coord['x'])
        pos_coords = np.array(emfisis_vars_coord['y'])
        if(np.nanmean(np.abs(pos_coords)) > 40):
            pos_coords = pos_coords / 6371.2
        x = pos_coords[:,0]
        y = pos_coords[:,1]
        z = pos_coords[:,2]

        pos_gsm = np.column_stack(( interp1(pos_times, x, times,np.nan), interp1(pos_times, y, times,np.nan), interp1(pos_times, z, times,np.nan) ))
        inds = np.where(np.logical_not(np.isnan( pos_gsm[:,0] )))
        times = times[inds]
        pos_gsm = pos_gsm[inds]
    else:
        times = np.array([])
        pos_gsm = np.array([[], [], []]).T
        inds = np.array([])
    return times, pos_gsm, inds


def get_fluxes_integrated_onR(mageis_vars, pos_coords, coord_onR, inds_notnan):
    times   = np.array(mageis_vars['FEDU']['x'])
    fluxes0 = np.array(mageis_vars['FEDU']['y'])
    alpha0  = np.array(mageis_vars['FEDU']['v1'])
    energy  = np.array(mageis_vars['FEDU']['v2'])
    times   = times[inds_notnan]
    fluxes0 = fluxes0[inds_notnan]
    fluxes0[np.where(fluxes0<0)] = 0
# #смотрим где последний бин по питч углу не поределен, и копируем туда данные из первого бина
#     inds = np.where(np.isnan(fluxes0[:,-1,:]))
#     fluxes0[(inds[0],-1*np.ones_like(inds[0]),inds[1])] = fluxes0[(inds[0],np.zeros_like(inds[0]),inds[1])]
# #тоже самое для -2
#     inds = np.where(np.isnan(fluxes0[:,-2,:]))
#     fluxes0[(inds[0],-2*np.ones_like(inds[0]),inds[1])] = fluxes0[(inds[0],np.ones_like(inds[0]),inds[1])]
#делаем сетку по alpha
    alpha = np.linspace(0, 180, 181)
    alpha = alpha[:-1] + np.diff(alpha)/2
#добавляем нулевой и 180 бин для интерполяции
#test chnge
    alpha01 = np.concatenate(([0], alpha0, [180]))
    fluxes01 = np.concatenate((fluxes0[:,[0],:], fluxes0, fluxes0[:,[-1],:]), axis=1)
    f = interp1d(alpha01, fluxes01, kind='linear', axis=1)
    fluxes = f(alpha) * np.mean(np.diff(alpha)) / np.mean(np.diff(alpha0))

    Bx0,By0,B0z = gp.ModelField(*pos_coords.T,*TFm.times_to_dateUt(times),Model='T96',CoordIn='GSM',CoordOut='GSM')
    Bx1,By1,Bz1 = gp.ModelField(*coord_onR.T,*TFm.times_to_dateUt(times),Model='T96',CoordIn='GSM',CoordOut='GSM')
    B0 = np.sqrt(Bx0**2 + By0**2 + B0z**2)
    B1 = np.sqrt(Bx1**2 + By1**2 + Bz1**2)

    alphaS = np.arcsin(np.outer(np.sqrt(B1 / B0), np.sin(alpha / 180 * np.pi))) * 180 / np.pi
    alphaS[:,np.where(alpha>90)] = 180 - alphaS[:,np.where(alpha>90)]

    fluxes[np.where(np.isnan(alphaS))] = np.nan
    fluxes = np.nansum(fluxes,axis=1)
    return fluxes #, alpha, energy


def get_fluxes_integrated_onRBSP(mageis_vars, inds_notnan):
    fluxes = np.array(mageis_vars['FEDU']['y'])
    fluxes = fluxes[inds_notnan]
    fluxes[np.where(fluxes<0)] = 0
    fluxes = np.nansum(fluxes,axis=1)
    return fluxes


def get_params(out_time, date):
    omni_time, dataOmni = omni.get_data(date)
    kp_time, kp_data    = Kp.get_data()

    Pdyn = interp1(omni_time, dataOmni['FlowPressure'], out_time)
    SymH = interp1(omni_time, dataOmni['SymH'],         out_time)
    By   = interp1(omni_time, dataOmni['ByGSM'],        out_time)
    Bz   = interp1(omni_time, dataOmni['BzGSM'],        out_time)  
    kp   = interp1(kp_time, kp_data, out_time)
    return Pdyn, SymH, By, Bz, kp


def get_data_flux_oneDay(date, R):
    """ 
    Get data projected at R for one day

    R: double -- radius of a sphere in RE at which projection should be

    date: datetime.date  --  date for projection 
    """
    dates_exclude = (datetime.date(2016, 10, 13),)
    if date in dates_exclude:
        return {}

    dirname = database+'data_flux/'
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    filename = dirname + 'data_flux_' + str(R) + '_' + date.strftime('%Y-%m-%d.dat')
    
    if not os.path.isfile(filename):
        date_current_str  = date.strftime('%Y-%m-%d/%H:%M:%S')
        date_next_str     = (date+datetime.timedelta(days=1)).strftime('%Y-%m-%d/%H:%M:%S')  
        mageis_vars = pyspedas.rbsp.mageis([date_current_str, date_next_str], probe='a', level='l3', rel='rel04', notplot=True)
        times = np.array(mageis_vars['FEDU']['x'])

        times, pos_coords, inds_notnan = get_coord_GSM(date_current_str, date_next_str, times)
        Trace = TFm.TraceFieldMy(date,times,pos_coords,Verbose=True)
        coord_onR1, coord_onR2 = get_coord_intersect(times, R, Trace)
        
        out_data = {}
        
        if(times.size):
            Pdyn, SymH, By, Bz, kp = get_params(times, date)
            out_data.update({'times':   times})
            out_data.update({'MLAT1':   interp1(Trace.time, Trace.MlatN, times)})
            out_data.update({'MLAT2':   interp1(Trace.time, Trace.MlatS, times)}) 
            out_data.update({'Lshell':  interp1(Trace.time, Trace.Lshell, times)})
            out_data.update({'fluxes1': get_fluxes_integrated_onR(mageis_vars, pos_coords, coord_onR1, inds_notnan)})
            out_data.update({'fluxes2': get_fluxes_integrated_onR(mageis_vars, pos_coords, coord_onR2, inds_notnan)})  
            out_data.update({'fluxesRBSP': get_fluxes_integrated_onRBSP(mageis_vars, inds_notnan)})
            out_data.update({'Pdyn':    Pdyn})
            out_data.update({'SymH':    SymH})
            out_data.update({'By':      By})
            out_data.update({'Bz':      Bz})
            out_data.update({'kp':      kp})

            # f1 = out_data.get('fluxes1')
            # f2 = out_data.get('fluxes2')
            # ft = out_data.get('fluxesRBSP')
            # f1e = np.sum(f1,axis=1)
            # f2e = np.sum(f2,axis=1)
            # fte = np.sum(ft,axis=1)
            # print(np.nanmax(f1e/fte))
            # print(np.nanmax(f2e/fte))
        
        with open(filename, 'wb') as fileData:
            pickle.dump(out_data,fileData,pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'rb') as fileData:
            out_data = pickle.load(fileData)
    return out_data


def calculate_dataflux(R, date_start_str="20120907", date_end_str="20191014",  days_step=30):
    """ 
    Get data projected at R from date_start to date_end

    R: double -- radius of a sphere in RE at which projection should be

    date_start: datetime.date  --  start date
    
    date_end: datetime.date  --  end date 

    days_step: uint  --  summing of dictionary routine interval in days
    """
    date_start  = datetime.date(int(date_start_str[0:4]), int(date_start_str[4:6]), int(date_start_str[6:8]))
    date_end    = datetime.date(int(date_end_str[0:4]), int(date_end_str[4:6]), int(date_end_str[6:8]))

    date_current = date_start
    one_day = datetime.timedelta(days=1)
    date_step = datetime.timedelta(days=days_step)
    
    while(date_current < date_end):
        date_end_p = min(date_current+date_step, date_end)
        while(date_current < date_end_p):
            get_data_flux_oneDay(date_current, R)
            date_current += one_day


def get_dataflux(R, date_start_str="20120907", date_end_str="20191014",  days_step=30):
    """ 
    Get data projected at R from date_start to date_end

    R: double -- radius of a sphere in RE at which projection should be

    date_start: datetime.date  --  start date
    
    date_end: datetime.date  --  end date 

    days_step: uint  --  summing of dictionary routine interval in days
    """
    date_start  = datetime.date(int(date_start_str[0:4]), int(date_start_str[4:6]), int(date_start_str[6:8]))
    date_end    = datetime.date(int(date_end_str[0:4]), int(date_end_str[4:6]), int(date_end_str[6:8]))

    date_current = date_start
    one_day = datetime.timedelta(days=1)
    date_step = datetime.timedelta(days=days_step)
    out_data = {}
    
    while(date_current < date_end):
        date_end_p = min(date_current+date_step, date_end)
        out_data_p = {}
        while(date_current < date_end_p):
            out_data_p = dict_conc(out_data_p, get_data_flux_oneDay(date_current, R))
            date_current += one_day
        out_data = dict_conc(out_data,out_data_p)

    return out_data


data = get_dataflux(4, "20121212", "20121213")

if __name__ == "__main__":
    import sys
    calculate_dataflux(int(sys.argv[1]), sys.argv[2], sys.argv[3])