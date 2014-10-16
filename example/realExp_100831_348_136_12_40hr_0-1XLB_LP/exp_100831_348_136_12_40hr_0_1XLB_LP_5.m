function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LP_5()
%% loading parameters
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
load('100831_348_136_12_40hr_0_1XLB_LP_input.mat');
snapPath = 'temp.mat';
param.CLUSTER_NUM = 32;
param.HES_FLAG = 0;
param.CLUSTER_NAME='killdevil1024_1';
param.D_HIST_PATH = 'DHist_l_9.05_t_0.40_p_14.55_NO_LATE_UPDATE.mat';
param.W_HIST_PATH = 'WHist_l_9.05_t_0.40_p_14.55_NO_LATE_UPDATE.mat';
param.LATE_UPDATE_FLAG = 0;
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], pDTemplate, [], elambda+3*elambda, etheta, ephi+2*ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LP_5_res.mat', 'expRec' );

end