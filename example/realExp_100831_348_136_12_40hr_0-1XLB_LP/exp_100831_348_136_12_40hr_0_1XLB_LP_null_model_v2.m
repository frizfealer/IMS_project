function [] = exp_100831_348_136_12_40hr_0_1XLB_LP_null_model_v2()
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/glmnet_matlab' ) );
addpath( genpath( '~/Ipopt_3118' ) );
addpath( genpath( '~/SLEP_package_4.1' ) );
addpath( genpath( '~/minFunc_2012' ) );
%% loading variables from exp 1
load( '100831_348_136_12_40hr_0_1XLB_LP_input_20150108.mat' );
[ param ] = setDLParameters();
lambda = 1e-32; theta = 1e-32; phi = 1e-32;
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], pDTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], 'temp.mat', param );
save( 'null_model_20150111_lf_iden.mat', expRec );
%%

