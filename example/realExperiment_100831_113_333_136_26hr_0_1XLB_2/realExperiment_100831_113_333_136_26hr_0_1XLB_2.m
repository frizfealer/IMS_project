function [ output_args ] = realExperiment_100831_113_333_136_26hr_0_1XLB_2(  tDataCubeFilePath, packageFolderPath, posIonTableFilePath, gridVec )
% %experimentTestingHParamsLimits The template of experiments on real data to
% %test the limits of hyper-parameters
% % tDataCubeFilePath = '100831_113_333_136_26hr_0-1XLB_1_inputs_env.mat';
% % packageFolderPath = your package path
% % posIonTableFilePath='molecule_profile_pos.csv';
% load( tDataCubeFilePath );
% 
% %% choose a tLen*tLen grid for testing hyper parameters
% tLen = 15;
% hLen = floor(tLen/2);
% blkSize = zeros(BlkDS.blkNum, 1);
% for i = 1:BlkDS.blkNum
%     blkSize(i) = length( BlkDS.B2GMap{i});
% end
% [~, cBlk] = max( blkSize );
% [i, j] = ind2sub( [size( tDataCube, 2 ) size( tDataCube, 3)], BlkDS.B2GMap{cBlk} );
% if median(i) - hLen  >= 1
%     upperY = median(i) - hLen;
% else
%     upperY = 1;
% end
% if median(i) + hLen <= max(i)
%     lowerY = median(i) + hLen;
% else
%     lowerY = max(i);
% end
% if median(j) - hLen >= 1
%     upperX = median(j) - hLen;
% else
%     upperX = 1;
% end
% if median(j) + hLen <= max(j)
%     lowerX = median(j) + hLen;
% else
%     lowerX = max(j);
% end
% ttDataCube = tDataCube( :, upperY:lowerY, upperX:lowerX );
% [ BlkDS ] = conBLKDS( ttDataCube );
% %% generate aMatrix, alpha matrix for leave-out testing.
% %do not leave out data, use all data
% aMatrix = ones( size( ttDataCube, 2 ), size( ttDataCube, 3) );
% %% Dictionary learning parameters setting
% snapPath = 'temp.mat';
% %setting up parameters
% param = [];
% param.OUTER_IT_NUM = 60; %number of  outer loop
% param.ADMM_IT_NUM = 100; %number of ADMM loop
% param.UP_D_IT_NUM = 200; %number of updating D loop
% param.HES_FLAG = 0; %whether using real Hessian when updating D
% param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
% param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
% param.LATE_UPDATE_FLAG = 1; %whether using late update on dictionary elements
% param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
% param.CLUSTER_NAME = 'killDevil1024'; %usually this is the default name
% param.CLUSTER_NUM = 32; %number of cluster uses, usually 12
% param.INIT_D = 'NNMF'; %the method of dictionary initialization
% 
% save( '100831_113_333_136_26hr_0_1XLB_2_inputs_env.mat', 'param', 'snapPath', 'BlkDS', 'aMatrix', 'tDataCube', ...
%     'ttDataCube', 'tDTemplate', 'tSpeciesM', 'tDIonName', 'gridVec', 'tMZAxis', '-v7.3' );

%% run on different hyper-parameters limits
load( '100831_113_333_136_26hr_0_1XLB_2_inputs_env.mat' );
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
param.CLUSTER_NAME = 'killdevil1024'; %usually this is the default name
param.CLUSTER_NUM = 32; %number of cluster uses, usually 12
gridLen = length( gridVec.fVec(:) );
expRecVec = cell(gridLen, 1);
for i = 1:gridLen
    cLambda = gridVec.fVec(i);
    cTheta = gridVec.sVec(i);
    cPhi = gridVec.tVec(i);
    [ expRec ] = dictionaryLearning_ADMM_v4( ttDataCube, [], tDTemplate, [], cLambda, cTheta, cPhi, aMatrix, BlkDS, [], snapPath, param );
    [ BICVal ] = computeBIC( 1, ttDataCube(:, :), expRec.outD, expRec.outW, expRec.outW0 );
    expRec.BICVal = BICVal;
    expRecVec{i} = expRec;
    save( 'realExperiment_100831_113_333_136_26hr_0_1XLB_2_result1.mat', 'expRecVec', '-v7.3' ) ;
end
BICsmall = inf;
gridsmall = 0;
for i = 1:gridLen
    if expRecVec{i}.BICVal < BICsmall
        BICsmall = expRecVec{i}.BICVal;
        gridsmall = i;
    end
end
cLambda = gridVec.fVec(gridsmall);
cTheta = gridVec.sVec(gridsmall);
cPhi = gridVec.tVec(gridsmall);

%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
[SLEN, IHEIGHT, IWIDTH] = size( tDataCube );
[ BlkDS ] = conBLKDS( tDataCube );
aMatrix = ones( size( tDataCube, 2 ), size( tDataCube, 3) );
param.OUTER_IT_NUM = 100; %number of  outer loop
[ expRecReal ] = dictionaryLearning_ADMM_v4( tDataCube, [], tDTemplate, [], cLambda, cTheta, cPhi, aMatrix, BlkDS, [], snapPath, param );
save( 'realExperiment_100831_113_333_136_26hr_0_1XLB_2_result2.mat', 'expRecReal', '-v7.3' ) ;
