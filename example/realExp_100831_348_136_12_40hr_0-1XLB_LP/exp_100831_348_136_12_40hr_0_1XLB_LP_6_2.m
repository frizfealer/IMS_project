function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LP_6_2()
%% loading parameters
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
load('100831_348_136_12_40hr_0_1XLB_LP_input.mat');
load('exp_100831_348_136_12_40hr_0_1XLB_LP_6_res.mat');
snapPath = 'temp_6.mat';
param.CLUSTER_NUM = 32;
param.HES_FLAG = 0;
param.OUTER_IT_NUM = 18;
param.CLUSTER_NAME='killdevil32';
param.D_HIST_PATH = 'DHist_l_11.32_t_0.40_p_14.55_NO_LATE_UPDATE_2.mat';
param.W_HIST_PATH = 'WHist_l_11.32_t_0.40_p_14.55_NO_LATE_UPDATE_2.mat';
param.LATE_UPDATE_FLAG = 0;
[ expRec_2 ] = dictionaryLearning_ADMM_v5( dataCube, expRec, pDTemplate, [], elambda+4*elambda, etheta, ephi+2*ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LP_6_2_res.mat', 'expRec_2' );

end
