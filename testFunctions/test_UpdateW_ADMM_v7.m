%% testing updateW_ADMM_v7 on synthetic data set 3.

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

%% different parameter testing.
[ WAdmmClass ] = initWADMM( InputData, size(simData.gD, 2) );
% run with default initialization
[n1] = updateW_ADMM_v7( InputData, simData.gD, WAdmmClass, 1 );
%% testing the case where D changes, and what RHO_INIT works better.
mD = simData.gD; mD(1,2) = 0.7;
n1.RHO_INIT = 1e-2;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-4;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-4 works better
n1.RHO_INIT = 1e-4;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-6;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-6 works better
mD = simData.gD; mD(1,2) = 0.9;
n1.RHO_INIT = 1e-2;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-4;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-4 works better
n1.RHO_INIT = 1e-4;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-6;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-6 works better
mD = simData.gD; mD(1,2) = 0.95;
n1.RHO_INIT = 1e-2;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-4;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-4 works better
n1.RHO_INIT = 1e-4;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-6;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-6 works better
mD = simData.gD; mD(1,2) = 0.99;
n1.RHO_INIT = 1e-2;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-4;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-4 works better
n1.RHO_INIT = 1e-4;
[n1_1] = updateW_ADMM_v7( InputData, mD, n1, 1 );
n1.RHO_INIT = 1e-6;
[n1_2] = updateW_ADMM_v7( InputData, mD, n1, 1 );
figure; plot( [n1_1.LPAry(:,1) n1_2.LPAry(:,1)]);
figure; plot( [n1_1.LPAry(:,2) n1_2.LPAry(:,2)]);
% 1e-6 works better

% setting different parameters in AdmmModel
% like AdmmModel.lambda = 1;
WAdmmClass.lambda=1;
%run another updates
WAdmmClass.RHO_INIT=1e-4;
[n1] = updateW_ADMM_v7( InputData, simData.gD, WAdmmClass, 1 );
WAdmmClass.RHO_INIT=1e-6;
[n2] = updateW_ADMM_v7( InputData, simData.gD, WAdmmClass, 1 );
%comparing the log likelihood of them
figure; plot( [n1.LPAry(:,1) n2.LPAry(:,1)]);

