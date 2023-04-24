












# from datetime import date
# import numpy as np
# import PlasmaFlowForecast as pff

# R = 3
# date_start = date(2012, 9, 7)
# date_end = date(2019, 10, 14)
 
# paramName = 'SymH'
# paramN = 30
# paramBins_b = paramN
# fluxesName = 'fluxes1'
# fluxesN = 100
# fluxesBins_b = np.logspace(-3, 8, fluxesN+1)

# energy          = pff.get_energy_rbsp_mageis()
# energyBins_N    = np.argwhere(energy>0)
# energyBins_N    = np.reshape(energyBins_N,(len(energyBins_N),))
# energyBins_N    = np.delete(energyBins_N, np.where(energyBins_N==13))

# paramBins_b  = pff.fit_fluxes(R, date_start, date_end, paramName, paramBins_b,
#                  fluxesName, fluxesBins_b, energyBins_N, is_save_figures=True, is_save_cvs=False)

# np.savetxt('fluxMean_log.cvs', paramBins_b, fmt='%.8e', delimiter=',')
# np.savetxt('fluxStd_log.cvs', paramBins_b, fmt='%.8e', delimiter=',')
# np.savetxt('paramBins_b.cvs', paramBins_b, fmt='%.8e', delimiter=',')

