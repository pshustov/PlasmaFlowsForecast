a = py.worktester.matlabfun();

fluxMean_log = double(a{1});
fluxStd_log  = double(a{2});
paramBins_b  = double(a{3});
energies     = double(a{4});

paramBins = paramBins_b(1:end-1) + diff(paramBins_b)/2;
paramBins_M = repmat(paramBins, 20, 1);
energies_M = repmat(energies', 1, 30);

errorbar(paramBins_M', fluxMean_log', fluxStd_log')
errorbar(paramBins_M(:,2:end-1)', fluxMean_log(:,2:end-1)', fluxStd_log(:,2:end-1)')


errorbar(energies_M, fluxMean_log, fluxStd_log)
set(gca,'XScale','log')


energies_log = log10(energies);
energies_log_M = log10(energies_M);

surf(paramBins_M, energies_log_M, fluxMean_log)



paramBins_M0 = paramBins_M(:,2:end-1);
energies_log_M0 = energies_log_M(:,2:end-1);
fluxMean_log0 = fluxMean_log(:,2:end-1);
fluxStd_log0 = fluxStd_log(:,2:end-1);

surf(paramBins_M0, energies_log_M0, fluxMean_log0)

np.log(gamma(np.arctan(a)/(np.pi/2)*k)*gamma(k - np.arctan(a)/(np.pi/2)*k)/gamma(k))