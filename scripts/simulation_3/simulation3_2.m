load( 'simulation3_env.mat' );
param.OUTER_IT_NUM = 100;
param.D_HIST_PATH = 'DHistCell_est1.mat';
param.W_HIST_PATH = 'WHistCell_est1.mat';
snapPath='simulation3_est1_temp.mat';
aMatrix = ones(30,30);
 [ sExpRec_est1 ] = dictionaryLearning_ADMM_v5( simData.gY, [], sDTemplate, [], 1e-6, 1e-6, 1e-6, aMatrix, BlkDS, [], snapPath, param );
 %[ mfitD, Matching ] = HungarianArrange( sExpRec_est1.outD, simData.gD(:, simData.usedElement) );
%figure; imagesc(mfitD(:,1:16)'*simData.gD(:, simData.usedElement)); colorbar;
[ elambda, ~, ephi ] = estimateHypParam( sExpRec_est1.outW, sExpRec_est1.outD, sDTemplate, BlkDS );
param.D_HIST_PATH = 'DHistCell_est2.mat';
param.W_HIST_PATH = 'WHistCell_est2.mat';
[ sExpRec_est2 ] = dictionaryLearning_ADMM_v5( simData.gY, [], sDTemplate, [], elambda, 1e-6, ephi, aMatrix, BlkDS, [], snapPath, param );
[ ~, etheta, ~ ] = estimateHypParam( sExpRec_est2.outW, sExpRec_est2.outD, sDTemplate, BlkDS );
param.D_HIST_PATH = 'DHistCell.mat';
param.W_HIST_PATH = 'WHistCell.mat';
param.OUTER_IT_NUM = 100;
snapPath='temp.mat';
[ sExpRec ] = dictionaryLearning_ADMM_v5( simData.gY, [], sDTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
save( 'simulation3RES.mat', 'sExpRec_est1', 'sExpRec_est2', 'sExpRec' );

[ mfitD, Matching ] = HungarianArrange( sExpRec.outD, simData.gD(:, simData.usedElement) );
gD = simData.gD(:, simData.usedElement);
for i = 1:16
    mfitD(:,i) = mfitD(:,i)/norm(mfitD(:,i));
    gD(:, i) = gD(:, i)/norm(gD(:, i));
end
r0 = mfitD(:,1:16)'*gD;
figure; 
imagesc(mfitD(:,1:16)'*gD); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 



%% compare with NNMF and pLSA
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
r1 = mfitD(:,1:16)'*gD;
figure;
imagesc(r1); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 


addpath D:\plsa2
AICVec = zeros(20, 1);
for i = 10:20
Par = []; Par.maxit = 200; Par.Leps = 1; Par.doplot = 0;
[Pw_z,Pd_z,Pz,Li] = pLSA_EM( simData.gY(:, BlkDS.indMap==1), i, Par );
for j = 1:i
    Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
end
D_pLSA = Pw_z;
[ AICVec(i) ] = computeAICc( simData.gY, D_pLSA, BlkDS, i);
end
%i = 16
[ mfitD, Matching ] = HungarianArrange( D_pLSA, simData.gD(:, simData.usedElement) );
r2 = mfitD(:,1:16)'*gD;
figure;
imagesc(r2); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 

figure;subplot(1,2,1);imagesc(r1);colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 
subplot(1,2,2);imagesc(r2); colorbar; ylabel( 'dictionary element from computation method' ); xlabel( 'dictionary element from ground truth' ); 

figure;subplot(2,2,1);imagesc(reshape( sum(simData.gY, 1), 30, 30 ) ); xlabel( 'width' ); ylabel('height');
figure;subplot(2,2,1);imagesc(gD'*gD); ylabel( 'element from ground truth' ); xlabel( 'element from ground truth' ); colorbar; 
subplot(2,2,2);imagesc(r0); colorbar; ylabel( 'element from DL' ); xlabel( 'element from ground truth' ); 
subplot(2,2,3);imagesc(r1);colorbar; ylabel( 'element from NNMF' ); xlabel( 'element from ground truth' ); 
subplot(2,2,4);imagesc(r2); colorbar; ylabel( 'element from pLSA' ); xlabel( 'element from ground truth' ); 

figure;subplot(1, 4,1);imagesc(reshape( sum(simData.gY, 1), 30, 30 ) ); xlabel( 'width' ); ylabel('height'); set(gca, 'square');
figure;ax = subplot(1, 4,1);imagesc(gD'*gD); ylabel( 'ground truth' ); xlabel( 'ground truth' ); colorbar;  axis square

figure; ha = tight_subplot(1,4,[.01 .03],[.1 .01],[.01 .01]); 
axes(ha(1));
imagesc(gD'*gD); ylabel( 'ground truth' ); xlabel( 'ground truth' ); colorbar;  axis square

subplot(1, 4,2);imagesc(r0); colorbar; ylabel( 'DL' ); xlabel( 'ground truth' );  axis square
subplot(1,4,3);imagesc(r1);colorbar; ylabel( 'NNMF' ); xlabel( 'ground truth' );   axis square
subplot(1,4,4);imagesc(r2); colorbar; ylabel( 'pLSA' ); xlabel( 'ground truth' );  axis square