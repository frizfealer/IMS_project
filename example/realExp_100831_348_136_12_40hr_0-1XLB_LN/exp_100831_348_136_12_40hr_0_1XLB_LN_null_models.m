function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LN_null_models()
%% loading variables from exp 1
load( '100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );
[ param ] = setDLParameters();
lambda = 1e-32; theta = 1e-32; phi = 1e-32;
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], nDTemplate2, [], lambda, theta, phi, aMatrix, BlkDS, [], 'temp.mat', param );
%%
save( 'null_model_20150111_lf_iden.mat', 'expRec' );
end