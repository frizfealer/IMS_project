%% testing dictionaryLearning_ADMM_v7 on synthetic data set 3.

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

%% setting the program
[ WA ] = initWADMM( InputData, size(simData.gD, 2) );
[ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
[ DL ] = initDL( WA, DI );
DL.SAVE_TEM_PERIOD = -1;
[ DL ] = dictionaryLearning_ADMM_v7( InputData, DL, 1,1, [] );