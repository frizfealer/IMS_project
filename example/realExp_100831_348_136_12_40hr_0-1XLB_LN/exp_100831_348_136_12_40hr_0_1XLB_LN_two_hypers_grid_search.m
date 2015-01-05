function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LN_two_hypers_grid_search( index )
%% loading variables from exp 1
load( '~/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );
% load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
[ param ] = setDLParameters();
[lVec, pVec] = meshgrid( lambdaVec , phiVec );
lambda = lVec(index);  phi = pVec(index);
tempName = ['exp_100831_348_136_12_40hr_0_1XLB_LN_tmp_' num2str(index) '_.mat' ];
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], nDTemplate, [], lambda, 1e-32, phi, aMatrix, BlkDS, [], tempName, param );
mkdir( '~/IMS_project-master/all_hypers_grid_search_results' );
resName = ['exp_100831_348_136_12_40hr_0_1XLB_LN_res_' num2str(index) '_.mat' ];
save( ['~/IMS_project-master/two_hypers_grid_search_results' resName], 'expRec');

