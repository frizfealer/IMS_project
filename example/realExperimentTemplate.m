function [ expRec ] = realExperimentTemplate( dataCubeFilePath, packageFolderPath, IonTableFilePathPos, IonTableFilePathNeg)
%experimentTemplate The template of experiments on real data
% dataCubeFilePath = 'dataCubeInfo.mat';
% packageFolderPath = your package path
% IonTableFilePath='molecule_profile_new_trim.csv';
load( dataCubeFilePath );
addpath( genpath( packageFolderPath ) );

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
intThreshold = 200;
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
%if aMatrix is set to ones e.g. aMatrix = ons( IHEIGHT, IWIDTH );,
%then there is no leave-out
%aMatrix: a logical matrix, as the same size of the grid, with 1 means the
%grids used in training process, and 0 means in testing process
[ aMatrix ] = genAMatrix( BlkDS, 0.2, 'byRow');


%% generate DTemplate
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
[ pDTemplate, pDIonName, pSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathPos, 1.007, 0.5 );
[ nDTemplate, nDIonName, nSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathNeg, 1.007, 0.5 );
save( 'exp1_inputs_env.mat' );

%% running dictionary learning
lambda = 1e-6; phi = 1e-6; theta = 1e-6;
snapPath = 'temp.mat';
%setting up parameters
param = [];
param.OUTER_IT_NUM = 2; %number of  outer loop
param.ADMM_IT_NUM = 1; %number of ADMM loop
param.UP_D_IT_NUM = 1; %number of updating D loop
param.HES_FLAG = 0; %whether using real Hessian when updating D
param.CLUSTER_NUM = 1; %number of cluster uses, usually 12
param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
param.LATE_UPDATE_FLAG = 1; %whether using late update on dictionary elements
param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
param.CLUSTER_NAME = 'local'; %usually this is the default name
param.INIT_D = 'NNMF'; %the method of dictionary initialization

%% need install the package SLEP and add the package to the path
% addpath( genpath( 'D:\SLEP_package_4.1\SLEP' ) );
[ expRec ] = dictionaryLearning_ADMM_v4( dataCube, [], DTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
[ elambda, etheta, ephi ] = estimateHypParam( expRec.outW, expRec.outD, DTemplate, BlkDS );
[ expRec ] = dictionaryLearning_ADMM_v4( dataCube, [], DTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
% expRec is a datastructure containing
% expRec.outD (dictionary)
% expRec.outW = outW (weight)
% expRec.outW0 = outW0; (weight offset)
% expRec.LPAry = LPAry; (Log posterior in outer loop)
% expRec.z0 = z0; expRec.z1 = z1; 
% expRec.u0 = u0; expRec.u1 = u1; (variables in ADMM)
% expRec.LPAryADMM = LPAryADMM; (Log posterior in ADMM loop)
% expRec.theta = theta; expRec.lambda = lambda; expRec.phi = phi; (hyperparameter in this experiment)
% expRec.rhoCell = rhoCell; expRec.resRecCell = resRecCell; (rho and residual history of the last ADMM loop)
save( 'expResult.mat', 'expRec' );

end

