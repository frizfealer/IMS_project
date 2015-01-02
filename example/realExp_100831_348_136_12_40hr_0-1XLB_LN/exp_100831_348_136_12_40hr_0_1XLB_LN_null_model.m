function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LN_null_model()
%% loading variables from exp 1
load( '100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );
[ param ] = setDLParameters();
lambda = 1e-32; theta = 1e-32; phi = 1e-32;
[ expRec ] = dictionaryLearning_ADMM_v6( dataCube, [], nDTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], 'temp.mat', param );
%%

end