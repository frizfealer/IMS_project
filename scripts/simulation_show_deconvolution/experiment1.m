%% deprecated, this scripts can not generate clear results
gD = eye(3); DTemplate = gD;
SLEN = 3; MLEN = 3; IHEIGHT = 20; IWIDTH = 20;
[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'diffusion', 1, [], 1, 1);
[ gY ] = genY_Poisson( gD, gW, gW0 );
BlkDS = conBLKDS( gY );

lambda = 1e-6; theta = 1e-6; phi = 1e-6;
aMatrix = ones( IHEIGHT, IWIDTH );
snapPath = 'temp.mat';
%setting up parameters
param = [];
param.OUTER_IT_NUM = 100; %number of  outer loop
param.ADMM_IT_NUM = 100; %number of ADMM loop
param.UP_D_IT_NUM = 200; %number of updating D loop
param.HES_FLAG = 0; %whether using real Hessian when updating D
param.CLUSTER_NUM = 12; %number of cluster uses, usually 12
param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
param.LATE_UPDATE_FLAG = 0; %whether using late update on dictionary elements
param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
param.LATE_UPDATE_INTTHRES = 0.8; %Under late update scenario, the percentage of intesity be explained used in the whole update process
param.CLUSTER_NAME = 'local'; %usually this is the default name
param.INIT_D = 'NNMF'; %the method of dictionary initialization
save( 'experiment1_set.mat' );

[ expRec ] = dictionaryLearning_ADMM_v4( gY, [], DTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
[ elambda, etheta, ephi ] = estimateHypParam( expRec.outW, expRec.outD, DTemplate, BlkDS );
[ expRec2 ] = dictionaryLearning_ADMM_v4( gY, [], DTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
param.INIT_D = 'random';
param.D_HIST_PATH = 'DHist_temp.mat';
[ expRec3 ] = dictionaryLearning_ADMM_v4( gY, [], DTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );