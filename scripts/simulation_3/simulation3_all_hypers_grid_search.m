function [ expRec, gridVec ] = simulation3_all_hypers_grid_search( index )
%% loading variables from exp 1
load( '/csbiohome01/ycharn/IMS_project-master/scripts/simulation_3/simulation3_20150108_lf_id.mat' );
% load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
[ param ] = setDLParameters();
param.OUTER_IT_NUM = 50;
[lVec, tVec, pVec] = meshgrid( lambdaVec...
    , thetaVec, phiVec );
lambda = lVec(index); theta = tVec(index); phi = pVec(index);
tempName = ['simulation3_tmp_' num2str(index) '_.mat' ];
fprintf( 'ready to run!\n' );
[ expRec ] = dictionaryLearning_ADMM_v6( simData.gY, [], sDTemplate, [], lambda, theta, phi, newAMatrix, BlkDS, [], tempName, param );
mkdir( '/lustre/scr/y/h/yharn/all_hypers_grid_search_results_simulation3' );
resName = ['simulation3_res_' num2str(index) '_.mat' ];
save( ['/lustre/scr/y/h/yharn/all_hypers_grid_search_results_simulation3/' resName], 'expRec');

