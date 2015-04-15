function [ ] = exp_100831_348_136_12_40hr_0_1XLB_LN_null_models_v2( varFilePath, nullModelPath )
%% loading variables from exp 1
load( varFilePath );
[ param ] = setDLParameters();
lambda = 1e-32; theta = 1e-32; phi = 1e-32;
aMatrix = ones( size( IMSD.BlkDS.indMap ) );
[ expRec ] = dictionaryLearning_ADMM_v6_2( IMSD.dataCube, [], DTemplate2, [], lambda, theta, phi, aMatrix, IMSD.BlkDS, [], 'temp.mat', param );
%%
save( nullModelPath, 'expRec' );
end