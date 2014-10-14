function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LN_4()
%% loading variables from exp 1
load('100831_348_136_12_40hr_0_1XLB_LN_input.mat');
snapPath = 'temp.mat';
param.CLUSTER_NUM = 12;
param.CLUSTER_NAME = 'local';
param.HES_FLAG = 0;
param.D_HIST_PATH = 'DHist_l_6.44_t_2.34_p_10.24.mat';
param.W_HIST_PATH = 'WHist_l_6.44_t_2.34_p_10.24.mat';
%% run an estimate of hyper parameters
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], nDTemplate, [], elambda, etheta, ephi + ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LN_4_res.mat', 'expRec' );
%%

end