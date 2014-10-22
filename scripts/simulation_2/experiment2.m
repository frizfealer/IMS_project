gD = zeros(3,3); gD(1, 1) = 2; gD(2, 1) = 1; gD(2, 2) = 2; gD(3, 2) = 1; gD(1, 3) = 1; gD(3, 3) = 2;
for i = 1:3
gD = gD ./ norm(gD(:,i));
end
SLEN = 3; MLEN = 3; IHEIGHT = 20; IWIDTH = 20;
%[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'diffusion', 1, [], 2, 1);
%% construct W specific for experiment 2
gW = zeros( MLEN, IHEIGHT, IWIDTH );
gW(1, (IHEIGHT-4+1):IHEIGHT, 1:4) = 1*abs(randn(4,4)+500);
midW = floor( IWIDTH / 2 );
gW(2, 1:4, (midW-2):(midW+1)) = 1*abs(randn(4,4)+300);
gW(3, (IHEIGHT-4+1):IHEIGHT, (IWIDTH-4+1):IWIDTH ) = 1*abs(randn(4,4)+500);
% 
% figure;
% for i = 1:3
% subplot(1,3,i); imagesc(reshape(gW(i,:),20,20));colorbar;
% end

h = fspecial( 'gaussian', 5 );
for i = 1:MLEN
    I = imfilter( reshape( gW(i,:), IHEIGHT, IWIDTH ), h );
    blurNum = 150;
    for j = 1:blurNum
        I = imfilter( reshape( I, IHEIGHT, IWIDTH ), h );
    end
    gW(i,:) = I(:);
end

figure;
for i = 1:3
subplot(1,3,i); imagesc(reshape(gW(i,:),20,20));colorbar;
end


[ gY ] = genY_Poisson( gD, gW, gW0 );
BlkDS = conBLKDS( gY );

for i = 1:3
    subplot(1, 3, i); imagesc( reshape( gY(i, :), 20, 20 ) ); colorbar;
end

lambda = 1e-6; theta = 1e-6; phi = 1e-6;
aMatrix = ones( IHEIGHT, IWIDTH );
snapPath = 'temp.mat';
%setting up parameters
param = setDLParameters();
param.OUTER_IT_NUM = 100; %number of  outer loop
param.ADMM_IT_NUM = 100; %number of ADMM loop
param.UP_D_IT_NUM = 200; %number of updating D loop
param.HES_FLAG = 0; %whether using real Hessian when updating D
param.CLUSTER_NUM = 4; %number of cluster uses, usually 12
param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
param.LATE_UPDATE_FLAG = 0; %whether using late update on dictionary elements
param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
param.LATE_UPDATE_INTTHRES = 0.8; %Under late update scenario, the percentage of intesity be explained used in the whole update process
param.CLUSTER_NAME = 'local'; %usually this is the default name
param.INIT_D = 'NNMF'; %the method of dictionary initialization
param.D_HIST_PATH = 'sDHist_expRec.mat';
param.W_HIST_PATH = 'sWHist_expRec.mat';

initVar = [];
DTemplate = zeros(3,3); DTemplate(1:2,1) = 1; DTemplate(2:3, 2) = 1; DTemplate([1,3],3) = 1;
initVar.outD = abs(randn(3,3));
initVar.outD(DTemplate==0) = 0;
for i = 1:3
    initVar.outD(:, i) = initVar.outD(:, i) / norm( initVar.outD(:, i) ); 
end
initVar.outW = zeros(3, 400);
initVar.outW0 = zeros( 20, 20 );
initVar.z1 = initVar.outW;
z0 = log(gY);
z0(z0==-inf)=0;
%if not in traing set, we shoud not initialize z0 according to it.
%take the average values
tmp = sum( z0(:, :), 2 ) / ( length( find( BlkDS.indMap == 1 ) ) );
for i = 1:400
    %if is in data area, but not in training set, 
    %set it to average spectrum
    if BlkDS.indMap(i) && ~aMatrix(i)
        z0(:, i) = tmp;
    end
end
initVar.z0 = z0(:, :);

figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape(gY(i,:),20,20) ); colorbar;
end

save( 'experiment2_set.mat' );


[ expRec0 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );

figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( expRec0.outW(i,:),20,20) ); colorbar;
end
[ elambda, ~, ephi ] = estimateHypParam( expRec0.outW, expRec0.outD, DTemplate, BlkDS );
param.OUTER_IT_NUM = 60;
[ expRec0_2 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], elambda, theta, ephi, aMatrix, BlkDS, [], snapPath, param );
[ ~, etheta, ~ ] = estimateHypParam( expRec0_2.outW, expRec0_2.outD, DTemplate, BlkDS );

[ expRec1 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
[ expRec3 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], lambda, etheta, phi, aMatrix, BlkDS, [], snapPath, param );
[ ~, etheta, ~ ] = estimateHypParam( expRec0.outW, expRec0.outD, DTemplate, BlkDS );
[ expRec4 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], lambda, etheta, phi, aMatrix, BlkDS, [], snapPath, param );
save('simulation_res.mat');

% expRec4 has OK result, and can be explained by our estimation method
%% compare with NNMF and pLSA
opt = statset('MaxIter',100,'Display','final');
[W, H] = nnmf( gY(:, BlkDS.indMap==1), 3, 'algorithm', 'als', 'options', opt, 'w0', initVar.outD, 'h0', initVar.outW(:, BlkDS.indMap==1) );
[W,H] = nnmf(gY(:, BlkDS.indMap==1), 3,'replicates',10,...
                   'options',opt,...
                   'algorithm','mult');
for i = 1:3
    W(:, i) = W(:, i) / norm( W(:, i) );
end
D_NNMF = W;
[ mfitD, Matching ] = HungarianArrange( D_NNMF, gD );

addpath D:\plsa2
AICVec = zeros(2, 1);
for i = 2:3
Par = []; Par.maxit = 200; Par.Leps = 1; Par.doplot = 0;
[Pw_z,Pd_z,Pz,Li] = pLSA_EM( gY(:, BlkDS.indMap==1), i, Par );
for j = 1:i
    Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
end
D_pLSA = Pw_z;
[ AICVec(i-1) ] = computeAICc( gY, D_pLSA, BlkDS, i);
end
[ mfitD, Matching ] = HungarianArrange( D_pLSA, gD(:, simData.usedElement) );
figure; imagesc(mfitD(:,1:16)'*gD(:, simData.usedElement)); colorbar;

%% draw bar plot
[sg, idx] = sort( gD(:), 'descend');
outMat = zeros( length(gD(:)), 4);
outMat(:,1) = sg;
[ mfitD, Matching ] = HungarianArrange( expRec4.outD, gD );
outMat(:,2) = mfitD(idx);
[ mfitD, Matching ] = HungarianArrange( D_NNMF, gD );
outMat(:,3) = mfitD(idx);
[ mfitD, Matching ] = HungarianArrange( D_pLSA, gD );
outMat(:,4) = [mfitD(idx(1:9))];
figure;bar(outMat)
legend('ground truth', 'DL', 'NNMF', 'pLSA');
xlabel('dictionary entries (sorted by intensity)');
ylabel('intensity');
export_fig figure3.pdf -transparent

