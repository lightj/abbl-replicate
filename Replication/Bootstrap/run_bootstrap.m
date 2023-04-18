function [] = run_bootstrap(arg1,arg2,workers,maxiterinc,maxitercons)

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

savepathincome = strcat('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Inc/ResultsNEWBOOT/inc_est',arg2,num2str(arg1));
savepathconsumption = strcat('/home/jdlight/ABBL_PMCMC/JOE_codes/Bootstrap/Cons/ResultsNEWBOOT/inc_est',arg2,num2str(arg1));

disp('GENERATING BOOTSTRAP SAMPLE')

% load data and generate hh_ids used in bootstrap sample
N = 2113;
data=load ('/home/jdlight/Data/data4est_vFINAL_WITH_DUAL_RESTRICTION.out');
household_ids = data(:,1);
household_ids = unique(household_ids);

N = size(household_ids,1);
bootstrap_household_ids = datasample(household_ids,N);

disp('RANDOM IDS:')
disp(bootstrap_household_ids)

% generate a new dataset using these household ids
boot_data = zeros(N*7,size(data,2));
for ii = 1:N
    ref_id = bootstrap_household_ids(ii);
    data_i = data(data(:,1)==ref_id,:);
    boot_data(1+(ii-1)*7:ii*7,:) = data_i;
end

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
