load('D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat');
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
[simData] = synthesizeData_Poisson( SLEN, MLEN, HEIGHT, WIDTH, sDTemplate, DOptions, WOptions, 1 );
BlkDS = conBLKDS(simData.gY);
save( 'simulation3_20150108.mat' );

%% run different algorithms
load( 'simulation3_20150108.mat' );
params = setDLParameters();
aMatrix = BlkDS.indMap;
snapPath = 'temp.mat';
 [ fExpRec ] = dictionaryLearning_ADMM_v6( simData.gY, [], sDTemplate, [], 1e-32, 1e-32, 1e-32, aMatrix, BlkDS, [], snapPath, params );
 %[ mfitD, Matching ] = HungarianArrange( sExpRec_est1.outD, simData.gD(:, simData.usedElement) );
%figure; imagesc(mfitD(:,1:16)'*simData.gD(:, simData.usedElement)); colorbar;
[ elambda, ~, ephi ] = estimateHypParam( fExpRec.W, fExpRec.D, sDTemplate, BlkDS, 0 );
% param.D_HIST_PATH = 'DHistCell_est2.mat';
% param.W_HIST_PATH = 'WHistCell_est2.mat';
 
[ sExpRec ] = dictionaryLearning_ADMM_v6( simData.gY, [], sDTemplate, [], elambda, 1e-32, ephi, aMatrix, BlkDS, [], snapPath, params );
[ ~, etheta, ~ ] = estimateHypParam( sExpRec.W, sExpRec.D, sDTemplate, BlkDS, 0 );

[ tExpRec ] = dictionaryLearning_ADMM_v6( simData.gY, [], sDTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, params );

save( 'simulation3_res_20150108.mat', 'fExpRec', 'sExpRec', 'tExpRec' );
[ gridVec ] = genHypParamGrid( [elambda etheta ephi], 6 );
aMatrix = leaveoutCBPatternsData( BlkDS, 10 );
wholeExp = [];
for i = 1:216
    [ expRec ] = dictionaryLearning_ADMM_v6( simData.gY, [], sDTemplate, [], gridVec.fVec(i), gridVec.sVec(i), gridVec.tVec(i), aMatrix, BlkDS, [], snapPath, params );
    wholeExp{i} = expRec;
    save( 'wholeExp.mat', 'wholeExp');
end


%% compare with our algorithm
[ mfitD, Matching ] = HungarianArrange( fExpRec.D, simData.gD );
gD = simData.gD;
for i = 1:MLEN
    mfitD(:,i) = mfitD(:,i)/norm(mfitD(:,i));
    gD(:, i) = gD(:, i)/norm(gD(:, i));
end
r0 = mfitD'*gD;
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
r0 = gD'*gD;
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
[ mfitD, Matching ] = HungarianArrange( D_NNMF, simData.gD(:, simData.usedElement) );
r1 = mfitD'*gD;
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
    Par = []; Par.maxit = 500; Par.Leps = 1; Par.doplot = 0;
    [Pw_z,Pd_z,Pz,Li] = pLSA_EM( simData.gY(:, BlkDS.indMap==1), i, Par );
    for j = 1:i
        Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
    end
    D_pLSA = Pw_z;
    [ AICVec(i) ] = computeAICc( simData.gY, Li(end), D_pLSA, BlkDS, i);
end
%i = 20
[ mfitD, Matching ] = HungarianArrange( D_pLSA, simData.gD );
r2 = mfitD'*gD;
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
