function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LP_all_hypers_grid_search_2( index )
%% loading variables from exp 1
load( '~/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LP/100831_348_136_12_40hr_0_1XLB_LP_input_20150108.mat' );
% load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
[ param ] = setDLParameters();
[lVec, tVec, pVec] = meshgrid( lambdaVec2...
    , thetaVec2, phiVec2 );
lambda = lVec(index); theta = tVec(index); phi = pVec(index);
tempName = ['exp_100831_348_136_12_40hr_0_1XLB_LP_2_tmp_' num2str(index) '_.mat' ];
fprintf( 'ready to run!\n' );
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], pDTemplate2, [], lambda, theta, phi, aMatrix, BlkDS, [], tempName, param );
mkdir( '/lustre/scr/y/h/yharn/all_hypers_grid_search_results_LP_2' );
resName = ['exp_100831_348_136_12_40hr_0_1XLB_LP_2_res_' num2str(index) '_.mat' ];
save( ['/lustre/scr/y/h/yharn/all_hypers_grid_search_results_LP_2/' resName], 'expRec');

