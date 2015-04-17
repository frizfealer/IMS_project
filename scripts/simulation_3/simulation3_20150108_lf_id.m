%% synthetic experiment 3: Dictionary recovery evaluation ISMB2015
% this experiment is shown in SUPPLEMENTARY DATA, and this script can generate 
% figure 1.
%% construct synthesized data for this experiment
% load('100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
% sDTemplate = nDTemplate; 
% sDTemplate = sDTemplate(:,1:20);
% tMZ = [];
% for i = 1:size(sDTemplate, 1)
%     if ~isempty( find( sDTemplate(i, :) ~= 0, 1 ) )
%         tMZ = [tMZ i];
%     end
% end
% sDTemplate = sDTemplate(tMZ, :);
% [SLEN, MLEN] = size(sDTemplate);
% HEIGHT = 30;
% WIDTH = 30;
% DOptions = [];
% DOptions.coheMax = 0.5;
% WOptions.type = 'diffusion';
% WOptions.supPrec = 1;
% WOptions.scale = 1e4;
% WOptions.sparsePrec = 0.8;
% [simData] = synthesizeData_Poisson( 'identity', SLEN, MLEN, HEIGHT, WIDTH, sDTemplate, DOptions, WOptions, 1 );
% BlkDS = conBLKDS(simData.gY);
% save( 'simulation3_20150108_lf_id.mat' );

%% run different algorithms
load( 'simulation3_20150108_lf_id.mat' );
params = setDLParameters();
aMatrix = BlkDS.indMap;
snapPath = 'temp.mat';
 [ fExpRec ] = dictionaryLearning_ADMM_v6( simData.gY, [], sDTemplate, [], 1e-32, 1e-32, 1e-32, aMatrix, BlkDS, [], snapPath, params );
maxLambda = ( max(simData.gY(:))*log(max(simData.gY(:))+1e-32)-max(simData.gY(:)) ) / max(simData.gY(:));
maxPhi =  ( max(simData.gY(:))*log(max(simData.gY(:))+1e-32)-max(simData.gY(:)) );
BlkDS = conBLKDS( simData.gY );
[sLen, hei, wid] = size( simData.gY );
% nLen = hei*wid;
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
maxTheta = ( max(simData.gY(:))*log(max(simData.gY(:))+1e-32)-max(simData.gY(:)) ) / sum(sum( abs( Rall*fExpRec.W' ) ));
lambdaVec = linspace( 1e-2, maxLambda, 5 )*1e-2;
phiVec = linspace( 1e-2, maxPhi, 5 )*1e-2;
thetaVec = linspace( 1e-2, maxTheta*2, 5 )*1e-2;
newAMatrix = leaveoutCBPatternsData( BlkDS, 10 );
save( 'simulation3_20150108_lf_id.mat' );

%% run with grid search
resFolder = 'all_hypers_grid_search_results_simulation3';
inputPath = 'simulation3_20150108_lf_id.mat';
[ lVec, tVec, pVec, valVec ] = readResultsAndTest( resFolder, inputPath );
tIdx = find( valVec == max(valVec) );
listing = dir( resFolder);
load( [resFolder '/' listing(tIdx+2)] );
bestExpRec = expRec;
sDTemplate2 = zeros( size( bestExpRec.D ) );
sDTemplate2(bestExpRec.D>1e-1)=1;

%% compare with our algorithm
[ mfitD, Matching ] = HungarianArrange( bestExpRec.D, simData.gD );
gD = simData.gD;
for i = 1:19
    mfitD(:,i) = mfitD(:,i)/norm(mfitD(:,i));
    gD(:, i) = gD(:, i)/norm(gD(:, i));
end
rMOLDL = mfitD(:,1:19)'*gD;

%% compare with gD
tD = simData.gD;
for i = 1:19
    tD(:, i) = tD(:, i) / norm(tD(:, i));
end
rgD = tD'*tD;

%% compare with NNMF
opt = statset('MaxIter',100,'Display','final');
% [W, H] = nnmf( simData.gY(:, BlkDS.indMap==1), 20, 'algorithm', 'als', 'options', opt );
[W,H] = nnmf(simData.gY(:, BlkDS.indMap==1), 20,'replicates',10,...
                   'options',opt,...
                   'algorithm','mult');
for i = 1:20
    W(:, i) = W(:, i) / norm( W(:, i) );
end
D_NNMF = W;
[ mfitD, Matching ] = HungarianArrange( D_NNMF, simData.gD );
rNNMF = mfitD(:,1:19)'*simData.gD;



%% compare with pLSA
AICVec = zeros(20, 1);
for i = 10:20
    Par = []; Par.maxit = 200; Par.Leps = 1; Par.doplot = 0;
    [Pw_z,Pd_z,Pz,Li] = pLSA_EM( simData.gY(:, BlkDS.indMap==1)+1e-32, i, Par );
    for j = 1:i
        Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
    end
    D_pLSA = Pw_z;
    [ AICVec(i) ] = computeAICc( simData.gY, Li(end), BlkDS, i);
end
i = 19;
Par = []; Par.maxit = 200; Par.Leps = 1; Par.doplot = 0;
[Pw_z,Pd_z,Pz,Li] = pLSA_EM( simData.gY(:, BlkDS.indMap==1)+1e-32, i, Par );

for j = 1:size(Pw_z, 2);
    Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
end
D_pLSA = Pw_z;
[ mfitD, Matching ] = HungarianArrange( D_pLSA, simData.gD );
rpLSA = mfitD'*simData.gD;

save( 'simulation3_20150108_lf_id_res.mat' );


%% draw figure 
% without running all algorithms, loading the experiment results and
% drawing the figure.
load( 'simulation3_20150108_lf_id_res.mat' );
figure; subplot(2, 2, 1);
[n, x] = hist( diag(rgD), 100 );
n = [1 1]; x = [1-0.005 1+0.005];
idx = find(~eye(size(rgD)));
[n2, x2] = hist( rgD(idx), 100 );
n2 = n2 / sum(n2);
bar(x, n, 'histc' );
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1.005]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );
h = title( 'A)' );
set(h, 'Units', 'normalized');
set(h, 'Position', [-0.02,1.05,0]);

subplot(2, 2, 2);
[n, x] = hist( diag(r1), 200 );
n = n / sum(n);
idx = find(~eye(size(rpLSA)));
[n2, x2] = hist( rpLSA(idx), round(length(rpLSA(idx))/2) );
n2 = n2 / sum(n2);
bar(x, n, 'histc');
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );
h = title( 'B)' );
set(h, 'Units', 'normalized');
set(h, 'Position', [-0.02,1.05,0]);
tmp = get( gca, 'position' ); tmp(1) = 0.53;
set( gca, 'position', tmp );

subplot(2, 2, 3);
[n, x] = hist( diag(rNNMF), 101 );
n = n / sum(n);
idx = find(~eye(size(rNNMF)));
[n2, x2] = hist( rNNMF(idx), 100 );
n2 = n2 / sum(n2);
bar(x, n, 'histc');
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1.005]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );
h = title( 'C)' );
set(h, 'Units', 'normalized');
set(h, 'Position', [-0.02,1.05,0]);

subplot(2, 2, 4);
[n, x] = hist( diag(rMOLDL), 100 );
n = n / sum(n);
idx = find(~eye(size(rMOLDL)));
[n2, x2] = hist( rMOLDL(idx), 100 );
n2 = n2 / sum(n2);
bar(x, n, 'histc');
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1.005]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );
h = title( 'D)' );
set(h, 'Units', 'normalized');
set(h, 'Position', [-0.02,1.05,0]);
tmp = get( gca, 'position' ); tmp(1) = 0.53;
set( gca, 'position', tmp );
