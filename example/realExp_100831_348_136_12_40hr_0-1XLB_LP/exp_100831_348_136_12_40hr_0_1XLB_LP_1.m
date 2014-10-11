function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LP_1()
%% loading parameters
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
load('100831_348_136_12,40hr_0-1XLB_LP.mzML_dc.mat');

%% trim data, only consdier m/z>=mzSmall and m/z <=mzLarge
mzSmall = 500;
mzLarge = max(mzAxis);
assert( mzSmall >= min( mzAxis ) );
assert( mzLarge <= max( mzAxis ) );
% dataCube, mzAxis is loaded from the *.mat file in dataCubeFilePath
tsIdx = find(mzAxis>=mzSmall, 1);
teIdx = find( mzAxis<=mzLarge , 1, 'last' );
dataCube = dataCube(tsIdx:teIdx, :, :);
mzAxis = mzAxis(tsIdx:teIdx);

%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
[SLEN, IHEIGHT, IWIDTH] = size( dataCube );
[ BlkDS ] = conBLKDS( dataCube );

%% trim data, only consider m/z channel has values larger then the value of  intThreshold accross all the grids
intThreshold = 100;
ins = dataCube(:,:);
ins = ins(:, BlkDS.indMap);
tmp = zeros( size( ins, 1 ), 1 );
for i = 1:length(tmp)
    tmp(i) = max( ins(i, :) );
end
tIdx = find(tmp>=intThreshold);
dataCube = dataCube(tIdx, :, :);
mzAxis = mzAxis(tIdx);

%% generate aMatrix, alpha matrix for leave-out testing.
%if aMatrix is set to ones e.g. aMatrix = ones( IHEIGHT, IWIDTH );,
%then there is no leave-out
%aMatrix: a logical matrix, as the same size of the grid, with 1 means the
%grids used in training process, and 0 means in testing process
aMatrix = ones( IHEIGHT, IWIDTH );

%% generate DTemplate
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
[ pDTemplate, pDIonName, pSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathPos, 1.007, 0.5 );

%% running dictionary learning
lambda = 1e-6; phi = 1e-6; theta = 1e-6;
snapPath = 'temp.mat';
%setting up parameters
[ param ] = setDLParameters();
param.CLUSTER_NAME = 'killdevil1024_1'; %usually this is the default name
param.CLUSTER_NUM = 64; %number of cluster uses, usually 12
param.D_ION_NAME = pDIonName;
param.D_HIST_PATH = 'DHist_l_0_t_0_p_0.mat';
param.W_HIST_PATH = 'WHist_l_0_t_0_p_0.mat';
save( '100831_348_136_12_40hr_0_1XLB_LP_input.mat' );

%% run an estimate of hyper parameters
[ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], pDTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
save( 'exp_100831_348_136_12_40hr_0_1XLB_LP_1_res.mat', 'expRec' );
[ elambda, etheta, ephi ] = estimateHypParam( expRec.outW, expRec.outD, pDTemplate, BlkDS );
clear expRec
save( '100831_348_136_12_40hr_0_1XLB_LP_input.mat' );
%%

end