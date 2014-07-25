function [ expRec, gridVec ] = realExperiment_100831_113_333_136_26hr_0-1XLB_1( posDataCubeFilePath, negDataCubeFilePath, packageFolderPath, posIonTableFilePath, negIonTableFilePath )
%experimentTemplate The template of experiments on real data
% posDataCubeFilePath = '100831_113_333_136_26hr_0-1XLB_LP.mat';
% negDataCubeFilePath = '100831_113_333_136_26hr_0-1XLB_LN.mat'
% packageFolderPath = 'D:\Users\YeuChern\GitHub\IMS_project';
% posIonTableFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_pos.csv';
% negIonTableFilePath = 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_neg.csv';
dataCubeVec = [];
mzAxisVec = [];
load( posDataCubeFilePath );
dataCubeVec{1} = dataCube;
mzAxisVec{1} = mzAxis;
load( negDataCubeFilePath );
dataCubeVec{2} = dataCube;
mzAxisVec{2} = mzAxis;
%addpath( genpath( packageFolderPath ) );
% addpath( genpath( 'D:\SLEP_package_4.1\SLEP' ) );

%% trim data, only consdier m/z>=mzSmall and m/z <=mzLarge
for i = 1:2
    dataCube = dataCubeVec{i};
    mzAxis = mzAxisVec{i};
    mzSmall = 500;
    mzLarge = max(mzAxis);
    assert( mzSmall >= min( mzAxis ) );
    assert( mzLarge <= max( mzAxis ) );
    % dataCube, mzAxis is loaded from the *.mat file in dataCubeFilePath
    tsIdx = find(mzAxis>=mzSmall, 1);
    teIdx = find( mzAxis<=mzLarge , 1, 'last' );
    dataCubeVec{i} = dataCube(tsIdx:teIdx, :, :);
    mzAxisVec{i} = mzAxis(tsIdx:teIdx);
end

%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
[~, h1, w1] = size( dataCubeVec{1} ); [~, h2, w2] = size( dataCubeVec{2} );
if h1~=h2 || w1 ~= w2
    fprintf( 'warning: the dimension of positve and negative data are not the sam.\n' );
end
[SLEN, IHEIGHT, IWIDTH] = size( dataCubeVec{1} );
[ BlkDS ] = conBLKDS( dataCubeVec{1} );

%% trim data, only consider m/z channel has values larger then the value of  intThreshold accross all the grids
for j = 1:2
    dataCube = dataCubeVec{j};
    mzAxis = mzAxisVec{j};
    intThreshold = 50;
    ins = dataCube(:,:);
    ins = ins(:, BlkDS.indMap);
    tmp = zeros( size( ins, 1 ), 1 );
    for i = 1:length(tmp)
        tmp(i) = max( ins(i, :) );
    end
    tIdx = find(tmp>=intThreshold);
    dataCubeVec{j} = dataCube(tIdx, :, :);
    mzAxisVec{j} = mzAxis(tIdx);
end

%% generate aMatrix, alpha matrix for leave-out testing.
%if aMatrix is set to ones e.g. aMatrix = ones( IHEIGHT, IWIDTH );,
%then there is no leave-out
%aMatrix: a logical matrix, as the same size of the grid, with 1 means the
%grids used in training process, and 0 means in testing process
aMatrix = ones( IHEIGHT, IWIDTH );

%% generate DTemplate
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
[ pDTemplate, pDIonName, pSpeciesM ] = genDTemplate( mzAxisVec{1}, posIonTableFilePath, 1.007, 0.5 );
[ nDTemplate, nDIonName, nSpeciesM ] = genDTemplate( mzAxisVec{2}, negIonTableFilePath, 1.007, 0.5 );
[ tSpeciesM, tDTemplate, tDIonName ] = ...
    mergePosNegTable( pSpeciesM, nSpeciesM, pDTemplate, nDTemplate, pDIonName, nDIonName, 0.5, negIonTableFilePath, mzAxisVec{2} );

%% combine positve and negative data cube
pCube = dataCubeVec{1};
nCube = dataCubeVec{2};
tDataCube = [pCube(:,:);nCube(:,:)];
SLEN = size( pCube, 1 ) + size( nCube, 1 );
tDataCube = reshape( tDataCube, SLEN, IHEIGHT, IWIDTH );
[ BlkDS ] = conBLKDS( tDataCube );
%% running dictionary learning
lambda = 1e-6; phi = 1e-6; theta = 1e-6;
snapPath = 'temp.mat';
%setting up parameters
param = [];
param.OUTER_IT_NUM = 30; %number of  outer loop
param.ADMM_IT_NUM = 100; %number of ADMM loop
param.UP_D_IT_NUM = 200; %number of updating D loop
param.HES_FLAG = 0; %whether using real Hessian when updating D
param.CLUSTER_NUM = 1; %number of cluster uses, usually 12
param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
param.LATE_UPDATE_FLAG = 1; %whether using late update on dictionary elements
param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
param.LATE_UPDATE_INTTHRES = 0.8; %Under late update scenario, the percentage of intesity be explained used in the whole update process
param.CLUSTER_NAME = 'local'; %usually this is the default name
param.INIT_D = 'NNMF'; %the method of dictionary initialization
save( '100831_113_333_136_26hr_0-1XLB_1_inputs_env.mat', 'param', 'snapPath', 'BlkDS', 'aMatrix', 'lambda', 'theta', 'phi', 'tDataCube', ...
    'tDTemplate', 'tSpeciesM', 'tDIonName', '-v7.3' );
%% need install the package SLEP and add the package to the path
% addpath( genpath( 'D:\SLEP_package_4.1\SLEP' ) );
[ expRec ] = dictionaryLearning_ADMM_v4( tDataCube, [], tDTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
[ elambda, etheta, ephi ] = estimateHypParam( expRec.outW, expRec.outD, DTemplate, BlkDS );
[ gridVec ] = genHypParamGrid( [elambda, etheta, ephi], 5 );
%elambda = 1.9303, etheta = 2.0707, ephi = 4.2295


end

