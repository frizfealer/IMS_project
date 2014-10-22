function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LP_9()
%% loading parameters
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
load('100831_348_136_12_40hr_0_1XLB_LP_input.mat');
snapPath = 'temp_9.mat';
param.CLUSTER_NUM = 32;
param.HES_FLAG = 0;
param.CLUSTER_NAME='killdevil32';
param.D_HIST_PATH = 'DHist_l_18.11_t_2.01_p_38.81_NO_LATE_UPDATE.mat';
param.W_HIST_PATH = 'WHist_l_18.11_t_2.01_p_38.81_NO_LATE_UPDATE.mat';
param.LATE_UPDATE_FLAG = 0;
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], pDTemplate, [], 8*elambda, 5*etheta, 8*ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LP_9_res.mat', 'expRec' );

end
