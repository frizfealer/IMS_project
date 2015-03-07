function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LP_basicSettings_v2()
%experiment scripts update @ 2015/03/05
%% load data, including dataCube and mzAxis
load( 'D:\IMS_DATA\100831_348_136_12,40hr_0-1XLB_LP\pp\100831_348_136_12,40hr_0-1XLB_LP.mzML_dc.mat' );
%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
IMSD = IMSData; clear IMSData;
IMSD.BlkDS = conBLKDS( IMSD.dataCube );

%% bin dataCube
[ mDataCube, mMZAxis, mappingFunc, ~, ~, ~ ] = binningDataCube( IMSD.dataCube, IMSD.mzAxis, IMSD.BlkDS, []);
fileID = fopen( 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_pos_v2.csv' );
C = textscan(fileID, '%f %s', 'delimiter', {','}); fclose(fileID); mpAry = C{1,1};
[ resMat, delta, resMat2, resMat3 ] = checkingBinFunction2( mappingFunc, IMSD.mzAxis, mpAry, mDataCube(:, IMSD.BlkDS.indMap==1),...
    IMSD.indMatrix, IMSD.dataCube(:, IMSD.BlkDS.indMap==1), [] );
ins = resMat3; ins(resMat>0.5) = 0; %figure; imagesc2(ins);
ins = max(ins, [], 2);
mBinIdx = ins >= 2;
[ IMSD.dataCube, mMZAxis, IMSD.mappingFunc, ~, ~, ~ ] = binningDataCube( IMSD.dataCube, IMSD.mzAxis, IMSD.BlkDS, mMZAxis(mBinIdx) );
IMSD.oMZAxis = IMSD.mzAxis; IMSD.mzAxis = mMZAxis; clear mMZAxis;
%checking
% [ resMat, delta, resMat2, resMat3 ] = checkingBinFunction2( mappingFunc, IMSD.mzAxis, mpAry, mDataCube(:, IMSD.BlkDS.indMap==1),...
%     IMSD.indMatrix, IMSD.dataCube(:, IMSD.BlkDS.indMap==1), [] );

%% retrieve data dimension
[SLEN, IHEIGHT, IWIDTH] = size( IMSD.dataCube );

%% trim data, only consider m/z channel has values larger then the value of  intThreshold accross all the grids
% intThreshold = 100;
% ins = IMSD.dataCube(:,:);
% ins = ins(:, IMSD.BlkDS.indMap);
% tmp = zeros( size( ins, 1 ), 1 );
% for i = 1:length(tmp)
%     tmp(i) = max( ins(i, :) );
% end
% hist(tmp, 100);
% tIdx = find(tmp>=intThreshold);
% rMZ = IMSD.mzAxis(tIdx);
% IMSD.dataCube = IMSD.dataCube(tIdx, :, :);
% IMSD.mzAxis = IMSD.mzAxis(tIdx);
% for i = 1:length(rMZ)
%     IMSD.mappingFunc(IMSD.mappingFunc == rMZ(i) ) = [];
%     IMSD.oMZAxis(IMSD.oMZAxis == rMZ(i) ) = [];
%     IMSD.indMatrix(IMSD.mappingFunc == rMZ(i), :) = [];
% end

%% generate aMatrix, alpha matrix for leave-out testing.
%if aMatrix is set to ones e.g. aMatrix = ones( IHEIGHT, IWIDTH );,
%then there is no leave-out
%aMatrix: a logical matrix, as the same size of the grid, with 1 means the
%grids used in training process, and 0 means in testing process
% aMatrix = ones( IHEIGHT, IWIDTH );
[ aMatrix ] = leaveoutCBPatternsData( BlkDS, 10 );
%% generate DTemplate
IonTableFilePathPos = 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_pos_v2.csv';
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
[ pDTemplate2, pDIonName2, pSpeciesM2 ] = genDTemplate_v4( IMSD.dataCube, IMSD.mzAxis, ...
    IonTableFilePathPos, 1.007, 0.5, IMSD.mappingFunc, IMSD.oMZAxis, IMSD.indMatrix);

save( '100831_348_136_12_40hr_0_1XLB_LP_input_20150306.mat' );

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


