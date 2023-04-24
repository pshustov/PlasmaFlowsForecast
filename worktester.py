import dataFlux

data1 = dataFlux.get_dataflux(4, "20140112", "20140113")
ft = data1.get('fluxesRBSP')
print(ft[1,:])
a= 1


# import datetime
# import numpy as np
# import PlasmaFlowForecast as pff
# from matplotlib import pyplot as plt

# def matlabfun():
#     R = 3
#     date_start = datetime.date(2012, 9, 7)
#     date_end = datetime.date(2019, 10, 14)
    
#     paramName = 'SymH'
#     paramN = 30
#     paramBins_b = paramN

#     fluxesName = 'fluxes1'
#     fluxesN = 100
#     fluxesBins_b = np.logspace(-3, 8, fluxesN+1)

#     energy = pff.get_energy_rbsp_mageis()
#     energyBins_N = np.argwhere(energy>0)
#     energyBins_N = np.reshape(energyBins_N,(len(energyBins_N),))
#     energyBins_N = np.delete(energyBins_N, np.where(energyBins_N==13))

#     fluxMean_log, fluxStd_log, paramBins_b  = pff.fit_fluxes_gauss(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N, is_save_figures=True)
#     paramBins = paramBins_b[:-1] + np.diff(paramBins_b) / 2

#     # fig = plt.figure(figsize=[6.5,19.5])
#     # ax = fig.subplots(1,1)
#     # for i in range(fluxMean_log.shape[0]):
#     #     ax.errorbar(paramBins[1:-1], fluxMean_log[i,1:-1], fluxStd_log[i,1:-1])
#     # plt.show()
#     return fluxMean_log, fluxStd_log, paramBins_b, energy[energyBins_N]
# # prt.get_equinumeros_subests_on_param('SymH', numOfBins=40, R=3)
# matlabfun()