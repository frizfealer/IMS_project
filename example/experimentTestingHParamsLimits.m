function [ output_args ] = experimentTestingHParamsLimits(  dataCubeFilePathPos, dataCubeFilePathNeg, packageFolderPath, IonTableFilePathPos, IonTableFilePathNeg )
%experimentTestingHParamsLimits The template of experiments on real data to
%test the limits of hyper-parameters
% dataCubeFilePath = 'dataCubeInfo.mat';
% packageFolderPath = your package path
% IonTableFilePath='molecule_profile_new_trim.csv';
if ~isempty( dataCubeFilePathPos )
    load( dataCubeFilePathPos );
end
if ~isempty( dataCubeFilePathNeg )
    load( dataCubeFilePathNeg );
end
addpath( genpath( packageFolderPath ) );
addpath( genpath( 'D:\SLEP_package_4.1\SLEP' ) );
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
%% choose a tLen*tLen grid for testing hyper parameters
tLen = 15;
hLen = floor(tLen/2);
blkSize = zeros(BlkDS.blkNum, 1);
for i = 1:BlkDS.blkNum
    blkSize(i) = length( BlkDS.B2GMap{i});
end
[~, cBlk] = max( blkSize );
[i, j] = ind2sub( [size( dataCube, 2 ) size( dataCube, 3)], BlkDS.B2GMap{cBlk} );
if median(i) - hLen  >= 1
    upperY = median(i) - hLen;
else
    upperY = 1;
end
if median(i) + hLen <= max(i)
    lowerY = median(i) + hLen;
else
    lowerY = max(i);
end
if median(j) - hLen >= 1
    upperX = median(j) - hLen;
else
    upperX = 1;
end
if median(j) + hLen <= max(j)
    lowerX = median(j) + hLen;
else
    lowerX = max(j);
end
tDataCube = dataCube( :, upperY:lowerY, upperX:lowerX );
[ BlkDS ] = conBLKDS( tDataCube );
%% generate aMatrix, alpha matrix for leave-out testing.
%do not leave out data, use all data
aMatrix = ones( size( tDataCube, 2 ), size( tDataCube, 3) );
%% generate DTemplate
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
modeCnt = 0;
if ~isempty(IonTableFilePathPos) 
    [ pDTemplate, pDIonName, pSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathPos, 1.007, 0.5 );
    modeCnt = modeCnt + 1;
end
if ~isempty(IonTableFilePathNeg)
    [ nDTemplate, nDIonName, nSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathNeg, 1.007, 0.5 );
    modeCnt = modeCnt + 3;
end
if modeCnt == 1
tDTemplate = pDTemplate;
end
if modeCnt == 3
    tDTemplate = nDTemplate;
end
if modeCnt == 4
     [ tSpeciesM, tDTemplate, tDIonName ] = mergePosNegTable( pSpeciesM, nSpeciesM, pDTemplate, nDTemplate, pDIonName, nDIonName, 0.5, IonTableFilePathNeg, inpeakMZ );
end
%% Dictionary learning parameters setting
snapPath = 'temp.mat';
%setting up parameters
param = [];
param.OUTER_IT_NUM = 60; %number of  outer loop
param.ADMM_IT_NUM = 100; %number of ADMM loop
param.UP_D_IT_NUM = 200; %number of updating D loop
param.HES_FLAG = 0; %whether using real Hessian when updating D
param.CLUSTER_NUM = 1; %number of cluster uses, usually 12
param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
param.LATE_UPDATE_FLAG = 1; %whether using late update on dictionary elements
param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
param.CLUSTER_NAME = 'local'; %usually this is the default name
param.INIT_D = 'NNMF'; %the method of dictionary initialization
%% run on different hyper-parameters limits
betaVec = 0.1:0.2:0.9;
lambdaAllVec = 1:0.5:5;
phiAllVec = 1:0.5:5;
maxLambdaAllVec = zeros( length( betaVec ), 1);
expRecVec = [];
for i=1:length(betaVec)
    cBeta = betaVec(i);
    for j = 1:length( lambdaAllVec )
        cLambdaAll = lambdaAllVec(j);
        lambda = cLambdaAll*cBeta;
        theta = cLambdaAll*(1-cBeta);
        [ expRec ] = dictionaryLearning_ADMM_v4( tDataCube, [], tDTemplate, [], lambda, theta, 0, aMatrix, BlkDS, [], snapPath, param );
        expRecVec = [expRecVec expRec];
        if max(expRec.outW(:)) < 1e-3
            maxLambdaAllVec(i) = lambdaAllVec(j);
        end
    end
end
for j = 1:length( lambdaAllVec )
    for i = 1:length( phiAllVec )
        cLambdaAll = lambdaAllVec(j);
        lambda = cLambdaAll*0.1;
        theta = cLambdaAll*0.9;
        cPhi = phiAllVec(i);
        [ expRec ] = dictionaryLearning_ADMM_v4( tDataCube, [], pDTemplate, [], lambda, theta, 0, aMatrix, BlkDS, [], snapPath, param );
    end
end

