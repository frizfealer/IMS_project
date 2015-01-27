function [ ] = null_models_LF_lG_L1()
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/glmnet_matlab' ) );
addpath( genpath( '~/Ipopt_3118' ) );
addpath( genpath( '~/SLEP_package_4.1' ) );
addpath( genpath( '~/minFunc_2012' ) );
%% loading variables from exp 1
load( '100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );
[ param ] = setDLParameters();
lambda = 1e-32; theta = 1e-32; phi = 1e-32;
param.D_CONSTRAINTS = 'L1';
param.linkFunc = 'log_gaussain';
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], DTemplate3, [], lambda, theta, phi, aMatrix, BlkDS, [], 'temp3.mat', param );
%%
save( 'null_model_20150111_lf_LG_DC_L1.mat', 'expRec' );
end