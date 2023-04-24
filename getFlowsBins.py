



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
    N1 = len(paramBins_b)-1
    N2 = len(fluxesBins_b)-1
    param = data.get(paramName)

    for eBin in energyBins_N:
        fluxes = data.get(fluxesName)[:,eBin]
        indsPositiveFluxes = np.where(fluxes>0)
        y = param[indsPositiveFluxes]
        fluxes = fluxes[indsPositiveFluxes]     
        for j in range(N2):   
            for i in range(N1):
                inds1 = np.logical_and(paramBins_b[i] < y, y <= paramBins_b[i+1])
                inds2 = np.logical_and(fluxesBins_b[j] < fluxes, fluxes <= fluxesBins_b[j+1])
                inds  = np.logical_and(inds1, inds2)
                bins[i,j] += np.sum(inds)
    return bins


def load_binsN(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N, N_load_step = 60, N_load_step_part = 20):
    """
    Function return bins array, with distribution on fluxes and "paramName", size(N_param, N_fluxes)

    R: double  --  radius of a projection sphere 

    date_start: datetime.date  --  start date of a projection data

    date_end: datetime.date  --  end date of a projection data

    paramName: str  --  name of parameter

    paramBins_b: ndaray  --  vector of boundaries for parameter bins 

    fluxesName: str  --  name of fluxes ("fluxes1" for north and "fluxes2" for south)

    fluxesBins_b: ndaray  --  vector of boundaries for fluxes bins

    energyBins_N: iterable  --  array of energy bins with int numbers
    """
    N1,N2    = len(paramBins_b)-1,len(fluxesBins_b)-1
    
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
        binsN    = np.zeros((N1,N2),dtype='int32')
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
