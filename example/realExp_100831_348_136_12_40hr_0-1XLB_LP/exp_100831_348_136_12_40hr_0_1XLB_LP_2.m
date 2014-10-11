function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LP_2()
%% loading parameters
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
load('100831_348_136_12,40hr_0-1XLB_LP_input.mat');
snapPath = 'temp.mat';
param.CLUSTER_NUM = 32;
param.HES_FLAG = 0;
param.D_HIST_PATH = 'DHist_l_2.26_t_0.40_p_4.85.mat';
param.W_HIST_PATH = 'WHist_l_2.26_t_0.40_p_4.85.mat';
%% run an estimate of hyper parameters
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], nDTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LP_2_res.mat', 'expRec' );

end