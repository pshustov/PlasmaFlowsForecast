import numpy as np
import pyspedas
import datetime
import pickle
import os
import hashlib
import json

from Globals import *
from PyGeopack.TraceField import TraceField
from math import ceil


def times_to_dateUt(times):
        pos_time_DT = [ datetime.datetime.utcfromtimestamp(times[i]) for i in range(len(times))]
        Date        = np.array([int(pos_time_DT[i].strftime('%Y%m%d')) for i in range(len(times))])
        ut          = np.array([(pos_time_DT[i] - datetime.datetime(pos_time_DT[i].year,pos_time_DT[i].month,pos_time_DT[i].day)).total_seconds() for i in range(len(times))]) / 3600
        return Date,ut


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


class TraceFieldMy(TraceField):
    def __init__(self, date, times = None, positions = None, Model='T96', CoordIn = 'GSM', CoordOut = 'GSM', 
                    alt=100.0, MaxLen=1000, DSMax=1.0,FlattenSingleTraces=True,Verbose=False,OutDtype='float64',**kwargs):
        
        t_step   = datetime.timedelta(minutes=2)
        allowed_names = ['time', 'pos', 'nstep', 'x', 'y', 'z', 
                    'Bx', 'By', 'Bz', 's', 'R', 'Rnorm', 
                    'GlatN','GlatS','MlatN','MlatS',
                    'GlonN','GlonS','MlonN','MlonS',
                    'GltN','GltS','MltN','MltS',
                    'Lshell','MltE','FlLen']


        filename  =self.__get_filename(date, times, positions, Model, CoordIn, CoordOut, 
                    alt, MaxLen, DSMax, FlattenSingleTraces, OutDtype, kwargs,
                    t_step, allowed_names)


        if not os.path.isfile(filename):
            if (times is None) and (positions is None):
                dateCurrentStr  = date.strftime('%Y-%m-%d/%H:%M:%S')
                dateNextStr     = (date+datetime.timedelta(days=1)).strftime('%Y-%m-%d/%H:%M:%S')
                emfisis_vars    = pyspedas.rbsp.emfisis([dateCurrentStr, dateNextStr], probe='a', level='l3', notplot=True, coord='gsm')

                times     = emfisis_vars.get('coordinates').get('x')
                positions = emfisis_vars.get('coordinates').get('y') / 6371.2 

            wherenan    = np.where(np.logical_not(np.isnan(positions[:,1])))
            times       = times[wherenan]
            positions   = positions[wherenan,:][0]

            pos_R       = np.sqrt(positions[:,0]**2 + positions[:,1]**2 + positions[:,2]**2)
            inds        = get_time_inds(times, t_step=t_step, R=pos_R)
            times       = times[inds]
            positions   = positions[inds,:]

            TraceField.__init__(self,*positions.T, *times_to_dateUt(times), Model, CoordIn, CoordOut, alt, MaxLen, DSMax,FlattenSingleTraces,Verbose,OutDtype,**kwargs)
            self.time   = times
            self.pos    = positions
            self.nstep  = self.nstep.astype('int32')

            if(self.__dict__):
                attr_names = list(self.__dict__.keys())
                for attr_name in attr_names:
                    if not any(attr_name == allowed_name for allowed_name in allowed_names):
                        exec('del self.%s' % (attr_name))

            with open(filename, 'wb') as fileTrace:
                pickle.dump(self.__dict__,fileTrace,pickle.HIGHEST_PROTOCOL)
        else:
            with open(filename, 'rb') as fileTrace:
                self.__dict__ = pickle.load(fileTrace)


    def __get_filename(self, date, times, positions, Model, CoordIn, CoordOut, alt, MaxLen, DSMax, FlattenSingleTraces, OutDtype, kwargs, t_step, allowed_names):
        # Verbose not indexing becous it just progress bar
        hashMD5 = hashlib.md5()    
        hashMD5.update(date.strftime('%Y%m%d').encode('utf-8'))
        hashMD5.update(np.array(times).astype('float16').tobytes())
        hashMD5.update(np.array(positions).astype('float16').tobytes())
        hashMD5.update(Model.encode('utf-8'))
        hashMD5.update(CoordIn.encode('utf-8'))
        hashMD5.update(CoordOut.encode('utf-8'))
        hashMD5.update(np.array(alt).astype('float16').tobytes())
        hashMD5.update(np.array(MaxLen).astype('int16').tobytes())
        hashMD5.update(np.array(DSMax).astype('float16').tobytes())
        hashMD5.update(np.array(FlattenSingleTraces).astype('bool_').tobytes())
        hashMD5.update(OutDtype.encode('utf-8'))
        hashMD5.update(json.dumps(kwargs).encode('utf-8'))
        hashMD5.update(np.array(t_step.total_seconds()).astype('float16').tobytes())
        hashMD5.update(''.join(allowed_names).encode('utf-8'))
        
        dirname = database+'rbsp_magtraces\\'
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        
        filename = dirname + 'magtrace_'+hashMD5.hexdigest()+'.dat'
        return  filename

