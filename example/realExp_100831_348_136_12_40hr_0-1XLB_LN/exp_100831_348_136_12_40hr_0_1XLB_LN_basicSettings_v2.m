function [] = exp_100831_348_136_12_40hr_0_1XLB_LN_basicSettings_v2( inputFilePath, outFilePath )
%experiment scripts, preparing the basic variables from xxx.mzML_dc.mat
%file
%% load data, including dataCube and mzAxis
%inputFilePath = 'D:\IMS_DATA\100831_348_136_12,40hr_0-1XLB_LN\pp\100831_348_136_12,40hr_0-1XLB_LN.mzML_dc.mat';
load( inputFilePath );

%% generate BlkDS, a data structure for bacteria community location information
IMSD = IMSData; clear IMSData;
IMSD.BlkDS = conBLKDS( IMSD.dataCube );

%% bin dataCube
[ mDataCube, mMZAxis, mappingFunc, ~, ~, ~ ] = binningDataCube( IMSD.dataCube, IMSD.mzAxis, IMSD.BlkDS, []);
IMSD.oMZAxis = IMSD.mzAxis; IMSD.mzAxis = mMZAxis; clear mMZAxis;
IMSD.dataCube = mDataCube; clear mDataCube;
IMSD.mappingFunc = mappingFunc; clear mappingFunc;
% fileID = fopen( 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_neg_v2.csv' );
% C = textscan(fileID, '%f %s', 'delimiter', {','}); fclose(fileID); mpAry = C{1,1};
% [ resMat, delta, resMat2, resMat3 ] = checkingBinFunction2( mappingFunc, IMSD.mzAxis, mpAry,...
%     IMSD.indMatrix, IMSD.dataCube(:, IMSD.BlkDS.indMap==1), [] );
% ins = resMat3; ins(resMat>0.5) = 0; %figure; imagesc2(ins);
% ins = max(ins, [], 2);
% mBinIdx = ins >= 2;
% [ IMSD.dataCube, mMZAxis, IMSD.mappingFunc, ~, ~, ~ ] = binningDataCube( IMSD.dataCube, IMSD.mzAxis, IMSD.BlkDS, mMZAxis(mBinIdx) );
%checking
% [ resMat, delta, resMat2, resMat3 ] = checkingBinFunction2( IMSD.mappingFunc, IMSD.oMZAxis, mpAry, mDataCube(:, IMSD.BlkDS.indMap==1),...
%     IMSD.indMatrix, IMSD.dataCube(:, IMSD.BlkDS.indMap==1), [] );

%% generate DTemplate
IonTableFilePathNeg = 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_neg_v2.csv';
%smallest molecule weight, set to H = 1.007
%m/z error +/- 0.5
% [ nDTemplate, nDIonName, nSpeciesM ] = genDTemplate( mzAxis, IonTableFilePathNeg, 1.007, 0.5 );
%[ DTemplate2, nDIonName2, nspeciesM2 ] = genDTemplate_v2( IMSD.dataCube, IMSD.mzAxis, IonTableFilePathNeg, 1.007, 0.5 );
% [ nDTemplate2, nDIonName2, nSpeciesM2 ] = genDTemplate_v4( IMSD.dataCube, IMSD.mzAxis, ...
%     IonTableFilePathNeg, 1.007, 0.5, IMSD.mappingFunc, IMSD.oMZAxis, IMSD.indMatrix);
[ DTemplate2, DIonName2, speciesM2 ] = genDTemplate_v3( IMSD.mzAxis, IonTableFilePathNeg, 5e-4 );
[ DTemplate2, DIonName2, SpeciesM2 ] = improveDTemplate( 'negative', DTemplate2, DIonName2, speciesM2 );
[ aMatrix ] = leaveoutCBPatternsData_v2( IMSD.BlkDS, 10 );
save( outFilePath );
% 
% [ maxLambda, minLambda ] = estimateLambdaMaxMin( IMSD.dataCube, expRec.D, expRec.W, expRec.W0, 'identity' );
% [ maxTheta, minTheta ] = estimateL2ThetaMaxMin( IMSD.dataCube, expRec.D, expRec.W, expRec.W0, 'identity' );
% [ maxPhi, minPhi ] = estimatePhiMaxMin( IMSD.dataCube, expRec.D, expRec.W, expRec.W0, 'identity' );
% lambdaVec = logspace( log10(minLambda), log10(maxLambda), 5 ); lambdaVec = [1e-32 lambdaVec];
% thetaVec = logspace( log10(minTheta), log10(maxTheta), 5 ); thetaVec = [1e-32 thetaVec];
% phiVec = logspace( log10(minPhi), log10(maxPhi), 5 ); phiVec = [1e-32 phiVec];
% 
% save( '100831_348_136_12_40hr_0_1XLB_LN_input_20150306.mat' );
% 
% [ aMatrix ] = leaveoutCBPatternsData_v2( IMSD.BlkDS, 10 );
% 
% save( '100831_348_136_12_40hr_0_1XLB_LN_input_20150306.mat' );

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