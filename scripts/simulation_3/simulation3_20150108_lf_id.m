load('/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
sDTemplate = nDTemplate; 
sDTemplate = sDTemplate(:,1:20);
tMZ = [];
for i = 1:size(sDTemplate, 1)
    if ~isempty( find( sDTemplate(i, :) ~= 0, 1 ) )
        tMZ = [tMZ i];
    end
end
sDTemplate = sDTemplate(tMZ, :);
[SLEN, MLEN] = size(sDTemplate);
HEIGHT = 30;
WIDTH = 30;
DOptions = [];
DOptions.coheMax = 0.5;
WOptions.type = 'diffusion';
WOptions.supPrec = 1;
WOptions.scale = 1e4;
WOptions.sparsePrec = 0.8;
[simData] = synthesizeData_Poisson( 'identity', SLEN, MLEN, HEIGHT, WIDTH, sDTemplate, DOptions, WOptions, 1 );
BlkDS = conBLKDS(simData.gY);
save( 'simulation3_20150108_lf_id.mat' );

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
resFolder = '';
inputPath = 'simulation3_20150108_lf_id.mat';
[ lVec, tVec, pVec, valVec ] = readResultsAndTest( resFolder, inputPath );
tIdx = find( valVec == max(valVec) );
listing = dir( resFolder);
load( [resFolder '/' listing(tIdx+2)] );
bestExpRec = expRec;
%actually index = 73 performs better?
sDTemplate2 = zeros( size( bestExpRec.D ) );
sDTemplate2(bestExpRec.D>1e-1)=1;


%% compare with our algorithm
[ mfitD, Matching ] = HungarianArrange( bestExpRec.D, simData.gD );
gD = simData.gD;
for i = 1:19
    mfitD(:,i) = mfitD(:,i)/norm(mfitD(:,i));
    gD(:, i) = gD(:, i)/norm(gD(:, i));
end
r0 = mfitD(:,1:19)'*gD;
% figure; 
%  set(gca,'fontsize',20, 'fontweight', 'bold')
% imagesc(r0); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 
[n, x] = hist( diag(r0), 100 );
n = n / sum(n);
idx = find(~eye(size(r0)));
[n2, x2] = hist( r0(idx), 100 );
n2 = n2 / sum(n2);
figure; bar(x, n, 'histc');
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1.005]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );
%% compare with gD
tD = simData.gD;
for i = 1:19
    tD(:, i) = tD(:, i) / norm(tD(:, i));
end
r0 = tD'*tD;
[n, x] = hist( diag(r0), 100 );
n = [1 1]; x = [1-0.005 1+0.005];
idx = find(~eye(size(r0)));
[n2, x2] = hist( r0(idx), 100 );
n2 = n2 / sum(n2);
figure; bar(x, n, 'histc' );
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1.005]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );



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
r1 = mfitD(:,1:19)'*simData.gD;
% figure;
% imagesc(r1); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 
[n, x] = hist( diag(r1), 101 );
n = n / sum(n);
idx = find(~eye(size(r1)));
[n2, x2] = hist( r1(idx), 100 );
n2 = n2 / sum(n2);
figure; bar(x, n, 'histc');
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1.005]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );


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
%i = 11
i = 19;
Par = []; Par.maxit = 200; Par.Leps = 1; Par.doplot = 0;
[Pw_z,Pd_z,Pz,Li] = pLSA_EM( simData.gY(:, BlkDS.indMap==1)+1e-32, i, Par );
for j = 1:i
    Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
end
D_pLSA = Pw_z;
[ mfitD, Matching ] = HungarianArrange( D_pLSA, simData.gD );
r2 = mfitD'*simData.gD;
% figure;
% imagesc(r2); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 
[n, x] = hist( diag(r1), 200 );
n = n / sum(n);
idx = find(~eye(size(r2)));
[n2, x2] = hist( r2(idx), round(length(r2(idx))/2) );
n2 = n2 / sum(n2);
figure; bar(x, n, 'histc');
hold on; h = bar(x2, n2, 'histc' );hold off;
set(h, 'facecolor', 'r');
legend( 'true dictionary elements', 'false dictionary elements' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'xlim', [0 1]);
xlabel( 'consie similarity', 'Fontsize', 20, 'Fontname','arial' );
ylabel( 'percentage', 'Fontsize', 20, 'Fontname','arial' );

save( 'simulation3_20150108_lf_id_res.mat' );

