%% testing updateW_ADMM_testing_v2 on synthetic data set 3.

%% preprocessing data
load( 'simulation3_20150108_lf_id' );
idx = find(BlkDS.indMap==1);
[I, J] = find(simData.gY(:, idx) == 0);
for i = 1:length(I)
	simData.gY(I(i), idx(J(i))) = 1;
end

InputData.dataCube =simData.gY;
InputData.BlkDS = BlkDS;
InputData.aMatrix = BlkDS.indMap;
[Rall] = genSparseGroupingMatrix( 30, 30, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end
InputData.Rblk = Rblk;
InputData.scaleFactor = 1;
InputData.logFY = logfactorial_e( simData.gY, 1e7);

%% run the ground-truth on all data.
[ WA ] = initWADMM( InputData, size(simData.gD, 2) );
[ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
[ gDL ] = initDL( WA, DI );
gDL.SAVE_TEM_PERIOD = -1;
[ gDL ] = dictionaryLearning_ADMM_v7( InputData, gDL, 0, 0, [] );


%% run the experiment on the one with held-out
[ aMatrix ] = leaveoutCBPatternsData_v2( BlkDS, 10 );
InputData.aMatrix = aMatrix;
[ WA ] = initWADMM( InputData, size(simData.gD, 2) );
[ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
[ hDL ] = initDL( WA, DI );
hDL.SAVE_TEM_PERIOD = -1;
[ hDL ] = dictionaryLearning_ADMM_v7( InputData, hDL, 0, 0, [] );
%run on the held-out subset
testData.dataCube= simData.gY(:, aMatrix==0 & BlkDS.indMap == 1 );
testData.BlkDS = conBLKDS( testData.dataCube );
testData.aMatrix = testData.BlkDS.indMap;
testData.scaleFactor = 1;
testData.logFY = logfactorial_e( testData.dataCube, 1e7 );
testData.Rblk = [];
[ TWA ] = initWADMM( testData, size(DL.D.D, 2) ); TWA.itNum = 200;
 TWA = updateW_ADMM_testing_v2( testData, DL.D.D, TWA, 1 );
 
 %% check if the estimated Ws in TWA and gDL are similar.
 targetIdx = find(aMatrix == 0 & BlkDS.indMap == 1);
 for i = 1:35
 figure; plot( gDL.W.W(:, targetIdx(i))); hold on; 
 plot( TWA.W(:,i), 'color', 'r' ); hold off;
 end
 
 %% testing choosing the best hyper parameters
% run the experiment on the one with held-out
load( 'simulation3_20150108_lf_id' );

[ InputData ] = initInputData( simData.gY, 10 );

% test W limits with lambda
lambdaVec = linspace(0.1, 10, 100);
lambdaVec = [0 lambdaVec];
AStruc = cell( length(lambdaVec), 1 );
parfor i = 1:length(lambdaVec)
    fprintf( '%d\n', i );
    [ WA ] = initWADMM( InputData, size(simData.gD, 2) );
    WA.lambda = lambdaVec(i);
    [AStruc{i}] = updateW_ADMM_v6( InputData, simData.uDTemplate, WA, 0 );
end
for i = 1:10:101
    ins = AStruc{i};
    insW = ins.W; insW = max(insW, [], 1 );
    figure;imagesc(reshape(insW, 30, 30 )); title(['lambda = ', num2str(ins.lambda)] );
end
% lambda should be 0, 1~5.

% test W limits with theta
thetaVec = logspace(-5, -2, 100);
thetaVec = [0 thetaVec];
AStruc = cell( length(thetaVec), 1 );
parfor i = 1:length(thetaVec)
    fprintf( '%d\n', i );
    [ WA ] = initWADMM( InputData, size(simData.gD, 2) );
    WA.theta = thetaVec(i);
    [AStruc{i}] = updateW_ADMM_v6( InputData, simData.uDTemplate, WA, 0 );
end
for i = 1:10:91
    ins = AStruc{i};
    insW = ins.W; insW = max(insW, [], 1 );
    figure;imagesc(reshape(insW, 30, 30 )); title(['theta = ', num2str(ins.lambda)] );
end
% theta should be 0, 2.0092e-05~6.5793e-04

% test D limts with phi
phiVec = logspace( 4, 8, 100 );
phiVec = [0 phiVec];
DStruc = cell( length(thetaVec), 1 );
parfor i = 1:length(thetaVec)
    fprintf( '%d\n', i );
    DI = initDInterior( simData.uDTemplate, InputData.dataCube,...
        [], [], [], [] );
    DI.phi = phiVec(i);
    [ DStruc{i}, ~ ] = updateD_v9_ipopt( InputData, ...
        AStruc{1}.W, AStruc{1}.W0, DI, 0 );
end
for i = 1:10:101
    ins = DStruc{i};
    insD = ins.D;
    figure;imagesc(insD); title(['phi = ', num2str(ins.phi)] );
end
% phi should be 0, 1e4~1e8

%run on the held-out subset
testCube= simData.gY(:, InputData.aMatrix==0 & InputData.BlkDS.indMap == 1 );
[ testData ] = initInputData( testCube, 0 );
[ TWA ] = initWADMM( testData, size(simData.uDTemplate, 2) ); TWA.itNum = 200;
lVec = linspace(1,5,5 ); lVec = [0 lVec];
tVec = logspace( log10(2.0092e-05), log10(6.5793e-04), 5 ); tVec = [0 tVec];
pVec = logspace( 4, 6, 5 ); pVec = [0 pVec];
[lv,tv,pv] = meshgrid( lVec, tVec, pVec );
lpVec = zeros( length(lv(:)), 1 );
parfor i = 1:length(lv(:))
    fprintf( '%d\n', i );
    [ WA ] = initWADMM( InputData, size(simData.gD, 2) );
    [ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
    WA.lambda = lv(i); WA.theta = tv(i); DI.phi = pv(i);
    [DL] = initDL( WA, DI );
    DL.SAVE_TEM_PERIOD = -1;
    [ nDL ] = dictionaryLearning_ADMM_v7( InputData, DL, 0, 0, [] );
    nTWA = updateW_ADMM_testing_v2( testData, nDL.D.D, TWA, 0 );
    lpVec(i) = nTWA.LPAry(find(nTWA.rhoAry~=0, 1, 'last'));
end
%load('lpvec.mat');
% i = 73 is the best within one standard deviation
i=73;
fprintf( '%d\n', i );
[ WA ] = initWADMM( InputData, size(simData.gD, 2) );
[ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
WA.lambda = lv(i); WA.theta = tv(i); DI.phi = pv(i);
[DL] = initDL( WA, DI );
DL.SAVE_TEM_PERIOD = -1;
[ nDL ] = dictionaryLearning_ADMM_v7( InputData, DL, 0, 0, [] );
nTWA = updateW_ADMM_testing_v2( testData, nDL.D.D, TWA, 0 );

