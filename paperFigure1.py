import numpy as np
import pyspedas
import PyGeopack as gp
import datetime
import TraceFieldMy as TFm
from scipy.interpolate import interp1d
from Globals import *

def interp1(x,y,xp):
    f = interp1d(x,y,bounds_error=False,fill_value=np.nan)
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
        x2[i] = interp1(trace_R[i,inds[i]:],trace_x[i,inds[i]:], R)
        y1[i] = interp1(trace_R[i,:inds[i]+1],trace_y[i,:inds[i]+1], R)
        y2[i] = interp1(trace_R[i,inds[i]:],trace_y[i,inds[i]:], R)
        z1[i] = interp1(trace_R[i,:inds[i]+1],trace_z[i,:inds[i]+1], R)
        z2[i] = interp1(trace_R[i,inds[i]:],trace_z[i,inds[i]:], R)

    return x1,y1,z1, x2,y2,z2


dirname = database+'data_paper\\figure1\\'

date = datetime.date(2016, 9, 7)
date_current_str  = date.strftime('%Y-%m-%d/%H:%M:%S')
date_next_str     = (date+datetime.timedelta(days=1)).strftime('%Y-%m-%d/%H:%M:%S') 
 
# data_rbsp = pyspedas.rbsp.mageis([date_current_str, date_next_str], probe='a', level='l3', rel='rel04', notplot=True)
data_rbsp = pyspedas.rbsp.emfisis([date_current_str, date_next_str], probe='a', level='l3', notplot=True, coord='gsm')
pos_time = data_rbsp.get('coordinates').get('x')
pos_position = data_rbsp.get('coordinates').get('y')
if(np.nanmean(np.abs(pos_position)) > 40):
    pos_position = pos_position / 6371.2

Trace = TFm.TraceFieldMy(date,pos_time,pos_position,Verbose=True)

trace_t = Trace.time
trace_x = Trace.x
trace_y = Trace.y
trace_z = Trace.z

R=4
X1,Y1,Z1, X2,Y2,Z2 = get_coord_intersect(pos_time, R, Trace)
R1 = np.column_stack((X1,Y1,Z1))
R2 = np.column_stack((X2,Y2,Z2))


np.savetxt(dirname+'times.csv',pos_time,fmt='%.8f',delimiter='\t')
np.savetxt(dirname+'position.csv',pos_position,fmt='%.8f',delimiter='\t')
np.savetxt(dirname+'trace_t.csv',trace_t,fmt='%.8f',delimiter='\t')
np.savetxt(dirname+'trace_x.csv',trace_x,fmt='%.8f',delimiter='\t')
np.savetxt(dirname+'trace_y.csv',trace_y,fmt='%.8f',delimiter='\t')
np.savetxt(dirname+'trace_z.csv',trace_z,fmt='%.8f',delimiter='\t') 
np.savetxt(dirname+'R1.csv',R1,fmt='%.8f',delimiter='\t') 
np.savetxt(dirname+'R2.csv',R2,fmt='%.8f',delimiter='\t') 