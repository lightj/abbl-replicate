function [] = run_par_bootstrap(arg1,arg2,workers,maxiterinc,maxitercons)

rng(int16(str2num(arg1)))

workers = int16(str2num(workers));
maxiterinc = str2num(maxiterinc)
maxitercons = str2num(maxitercons)

% set the parallel profile as local
pc = parcluster('local')
% set the temp file location to MATLABWORKDIR variable set in submit script
pc.JobStorageLocation = strcat(getenv('MATLABWORKDIR'))
% open a pool of 12 worker processes
disp('workers...')
disp(workers)
disp('type of workers...')
disp(class(workers))
disp('class of maxiterinc...')
disp(class(maxiterinc))
disp(maxitercons)
parpool(pc,workers)

savepathincome = strcat('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Inc_parametric/ResultsNEWBOOT/inc_est',arg2,num2str(arg1));
savepathconsumption = strcat('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Cons_parametric/ResultsNEWBOOT/inc_est',arg2,num2str(arg1));

% path that allows for the parametric sampler
addpath('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/')

% we use the same bootstrap estimation functions as the non-parametric bootstrap
disp('GENERATING BOOTSTRAP SAMPLE')

% simulate panel using estimated parameters
boot_data = simulate_panel();

% estimate income parameters
disp('ESTIMATING INCOME PARAMETERS')
SMC_n_vF(boot_data,maxiterinc,savepathincome)

% estimate consumption parameters
disp('ESTIMATING CONSUMPTION PARAMETERS')
PMCMC_simple_cons_vF(boot_data,maxitercons,savepathincome,savepathconsumption)

% close the workers and exit
poolobj = gcp('nocreate');
delete(poolobj);
end
