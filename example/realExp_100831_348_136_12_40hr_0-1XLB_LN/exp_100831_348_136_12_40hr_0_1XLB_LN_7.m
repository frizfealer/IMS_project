function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LN_7()
%% loading variables from exp 1
load('100831_348_136_12_40hr_0_1XLB_LN_input.mat');
snapPath = 'temp_7.mat';
param.CLUSTER_NUM = 12;
param.CLUSTER_NAME = 'local';
param.HES_FLAG = 0;
param.LATE_UPDATE_FLAG = 0;
param.D_HIST_PATH = 'DHist_l_5.81_t_7.68_p_6.58.mat';
param.W_HIST_PATH = 'WHist_l_5.81_t_7.68_p_6.58.mat';
%% using new set of estimated parameters
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], nDTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LN_7_res.mat', 'expRec' );
%%

end