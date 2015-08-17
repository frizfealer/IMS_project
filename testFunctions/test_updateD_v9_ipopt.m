%% testing updateD_v9_ipopt on synthetic data set 3.

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
[ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
DI.phi = 0;
[ DI, ~ ] = updateD_v9_ipopt( InputData, simData.gW, simData.gW0, ...
    DI, 1 );
figure; subplot(1,2,1); imagesc(simData.gD); subplot(1,2,2); imagesc(DI.D);

%test with ph = 100
[ DI ] = initDInterior( simData.uDTemplate, ...
    InputData.dataCube, [], [], [], [] );
DI.phi = 1e7;
[ DI, ~ ] = updateD_v9_ipopt( InputData, simData.gW, simData.gW0, ...
    DI, 1 );
figure; subplot(1,2,1); imagesc(simData.gD); subplot(1,2,2); imagesc(DI.D);
