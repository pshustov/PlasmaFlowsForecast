function matlabTest()
R = 3;

date_start = datetime(2012,9,7);
date_start_vec = int32(datevec(date_start));
date_start = py.datetime.date(date_start_vec(1), date_start_vec(2), date_start_vec(3));

date_end = datetime(2019,10,14);
date_end_vec = int32(datevec(date_end));
date_end = py.datetime.date(date_end_vec(1), date_end_vec(2), date_end_vec(3));

paramName = py.str('SymH');
paramN = int32(30);
paramBins_b = py.numpy.linspace(int32(-200),int32(50),paramN+1);

fluxesName = py.str('fluxes1');
fluxesN = int32(40);
fluxesBins_b = py.numpy.logspace(int32(0), int32(6), fluxesN+1);

energyBins_N = py.numpy.array(py.range(int32(5),int32(21)));

dataGauss_py = py.projectionAtR.fit_fluxes_gauss(R, date_start, date_end, paramName, paramBins_b, fluxesName, fluxesBins_b, energyBins_N);

fluxMean_log = double(dataGauss_py{1});
fluxStd_log = double(dataGauss_py{2});

paramBins_b = double(paramBins_b);
fluxesBins_b = double(fluxesBins_b);
energyBins_N = double(energyBins_N);

end

