function [ expRec ] = exp_100831_348_136_12_40hr_0_1XLB_LN_all_hypers_grid_search_v3(index)
% INPUT ARGUMENTS:
%   index, index in the hyperparameters grid searching.
projectPath = '~/IMS_project-master';
varPath = 'example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_vars.mat';
outFolderPath = '/lustre/scr/y/h/yharn/all_hypers_grid_search_results_neg';
resPath = 'exp_100831_348_136_12_40hr_0_1XLB_LN_res_v3_';
load( [projectPath '/' varPath] );
% load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
[ param ] = setDLParameters();
[lVec, tVec, pVec] = meshgrid( lambdaVec...
    , thetaVec, phiVec );
lambda = lVec(index); theta = tVec(index); phi = pVec(index);
tempName = ['exp_100831_348_136_12_40hr_0_1XLB_LN_3_tmp_' num2str(index) '_.mat' ];
fprintf( 'ready to run!\n' );
[ expRec ] = dictionaryLearning_ADMM_v6_2( IMSD.dataCube, [], DTemplate2, [], lambda, theta, phi, aMatrix, IMSD.BlkDS, [], tempName, param );
mkdir( outFolderPath );
resName = [ resPath num2str(index) '_.mat' ];
save( [outFolderPath '/' resName], 'expRec');

