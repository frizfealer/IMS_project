function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LP_4()
%% loading parameters
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
load('100831_348_136_12_40hr_0_1XLB_LP_input.mat');
snapPath = 'temp.mat';
param.CLUSTER_NUM = 32;
param.HES_FLAG = 0;
param.CLUSTER_NAME='killdevil1024_2';
param.D_HIST_PATH = 'DHist_l_3.39_t_0.40_p_7.27_NO_LATE_UPDATE.mat';
param.W_HIST_PATH = 'WHist_l_3.39_t_0.40_p_7.27_NO_LATE_UPDATE.mat';
param.LATE_UPDATE_FLAG = 0;
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], pDTemplate, [], elambda+elambda/2, etheta, ephi+ephi/2, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LP_4_res.mat', 'expRec' );

end