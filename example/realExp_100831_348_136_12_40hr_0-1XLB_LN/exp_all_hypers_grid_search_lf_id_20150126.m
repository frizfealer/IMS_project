function [ expRec, gridVec ] = exp_all_hypers_grid_search_lf_id_20150126( index )
%% loading variables from exp 1
load( '~/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );
% load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
[ param ] = setDLParameters();
param.linkFunc = 'identity';
[lVec, tVec] = meshgrid( lambdaVec_lf_id...
    , thetaVec_lf_id );
lambda = lVec(index); theta = tVec(index);
tempName = ['exp_100831_348_136_12_40hr_0_1XLB_LN_tmp_' num2str(index) '_.mat' ];
fprintf( 'ready to run!\n' );
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], DTemplate3, [], lambda, theta, 1e-32, aMatrix, BlkDS, [], tempName, param );
mkdir( '/lustre/scr/y/h/yharn/all_hypers_grid_search_results_lf_id_20150126' );
resName = ['exp_100831_348_136_12_40hr_0_1XLB_LN_res_' num2str(index) '_.mat' ];
save( ['/lustre/scr/y/h/yharn/all_hypers_grid_search_results_lf_id_20150126/' resName], 'expRec');

