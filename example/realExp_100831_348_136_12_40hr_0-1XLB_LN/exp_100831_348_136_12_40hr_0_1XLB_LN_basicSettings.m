function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LN_basicSettings()
%experiment scripts update @ 2014/12/31
%% load data, including dataCube and mzAxis
load( 'D:\IMS_DATA\100831_348_136_12,40hr_0-1XLB_LN\100831_348_136_12,40hr_0-1XLB_LN_noBin\100831_348_136_12,40hr_0-1XLB_LN.mzML_dc.mat' );
%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
BlkDS = conBLKDS( dataCube );
%% bin dataCube
[ mDataCube, mMZAxis, ~ ] = binningDataCube( dataCube, mzAxis, BlkDS );
dataCube = mDataCube; mzAxis = mMZAxis;
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

%% retrieve data dimension
[SLEN, IHEIGHT, IWIDTH] = size( dataCube );

%% trim data, only consider m/z channel has values larger then the value of  intThreshold accross all the grids
intThreshold = 10;
ins = dataCube(:,:);
ins = ins(:, BlkDS.indMap);
tmp = zeros( size( ins, 1 ), 1 );
for i = 1:length(tmp)
    tmp(i) = max( ins(i, :) );
end
hist(tmp, 100);
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
IonTableFilePathNeg = 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_neg.csv';
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
[ nDTemplate, nDIonName, nSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathNeg, 1.007, 0.5 );

save( '100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat' );

% %% running dictionary learning
% lambda = 1e-6; phi = 1e-6; theta = 1e-6;
% snapPath = 'temp_1.mat';
% %setting up parameters
% [ param ] = setDLParameters();
% param.D_ION_NAME = nDIonName;
% % param.D_HIST_PATH = 'DHist_l_0_t_0_p_0.mat';
% % param.W_HIST_PATH = 'WHist_l_0_t_0_p_0.mat';
% save( '100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
% 
% %% run an estimate of hyper parameters
% [ expRec ] = dictionaryLearning_ADMM_v5( dataCube, [], nDTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
% save( 'exp_100831_348_136_12_40hr_0_1XLB_LN_1_res.mat', 'expRec' );
% [ elambda, ~, ephi ] = estimateHypParam( expRec.outW, expRec.outD, nDTemplate, BlkDS );
% clear expRec
% [ expRec1_2 ] = dictionaryLearning_ADMM_v5( dataCube, [], nDTemplate, [], elambda, theta, ephi, aMatrix, BlkDS, [], snapPath, param );
% save( 'exp_100831_348_136_12_40hr_0_1XLB_LN_1_2_res.mat', 'expRec1_2' );
% [ ~, etheta, ~ ] = estimateHypParam( expRec1_2.outW, expRec1_2.outD, nDTemplate, BlkDS );
% clear expRec1_2
% save( '100831_348_136_12_40hr_0_1XLB_LN_input.mat' );
% %%

end