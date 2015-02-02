gD = zeros(3,3); gD(1, 1) = 2; gD(2, 1) = 1; gD(2, 2) = 2; gD(3, 2) = 1; gD(1, 3) = 1; gD(3, 3) = 2;
for i = 1:3
gD = gD ./ norm(gD(:,i), 1);
end
SLEN = 3; MLEN = 3; IHEIGHT = 20; IWIDTH = 20;
%[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'diffusion', 1, [], 2, 1);
%% construct W specific for this experiment
supPerc = 1;
sparsePerc = 1;
scale = [1000 1000 1000];
[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'diffusion', supPerc, sparsePerc, scale, 0);
gW0(gW(1,:)~=0)=0;
[ gY ] = genY_differentNoise( 'identity', gD, gW, gW0 );
BlkDS = conBLKDS( gY );

figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( gY(i, :), 20, 20 ) ); colorbar;
end

figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( gW(i, :), 20, 20 ) ); colorbar;
end

lambda = 1e-32; theta = 1e-32;
aMatrix = ones( IHEIGHT, IWIDTH );
snapPath = 'temp.mat';
%setting up parameters
param = setDLParameters();

initVar = [];
DTemplate = zeros(3,3); DTemplate(1:2,1) = 1; DTemplate(2:3, 2) = 1; DTemplate([1,3],3) = 1;
% initVar.outD = abs(randn(3,3));
% initVar.outD(DTemplate==0) = 0;
% for i = 1:3
%     initVar.outD(:, i) = initVar.outD(:, i) / norm( initVar.outD(:, i) ); 
% end
% initVar.outW = zeros(3, 400);
% initVar.outW0 = zeros( 20, 20 );
% initVar.z1 = initVar.outW;
% z0 = log(gY);
% z0(z0==-inf)=0;
% %if not in traing set, we shoud not initialize z0 according to it.
% %take the average values
% tmp = sum( z0(:, :), 2 ) / ( length( find( BlkDS.indMap == 1 ) ) );
% for i = 1:400
%     %if is in data area, but not in training set, 
%     %set it to average spectrum
%     if BlkDS.indMap(i) && ~aMatrix(i)
%         z0(:, i) = tmp;
%     end
% end
% initVar.z0 = z0(:, :);
% 
% figure;
% for i = 1:3
%     subplot(1, 3, i); imagesc( reshape(gY(i,:),20,20) ); colorbar;
% end
% 
save( 'experiment_fusion_env_lf_id.mat' );
%testing
iD = initD(gY, DTemplate, 'NNMF', [], 1 );
[ Dbest, ~ ]= updateD_v8_ipopt( 'identity', 'L1', gY, gW, gW0, iD, DTemplate, BlkDS.indMap, 0, 1e-32, 1, 200, 5e-3, [] );
%can get the correct results!
W_TOL=5e-3; D_LOWER_BOUND = 1e-2;
[WResStruct] = updateW_ADMM_v4( 'identity', gY, gD, aMatrix, 100, 0, 0, 0, [], [], 1, 0, W_TOL, D_LOWER_BOUND, [], [] );
for i = 1:3
    figure;
    subplot(1,3, 1); imagesc(reshape(WResStruct.W(i,:), 20, 20 ) ); colorbar; 
    subplot(1,3, 2); imagesc(reshape(gW(i,:), 20, 20 ) ); colorbar;
    subplot(1,3 ,3); imagesc(reshape(abs(WResStruct.W(i,:)-gW(i,:)), 20, 20 ) ); colorbar;
end
[WResStruct] = updateW_ADMM_v4( 'identity', gY, iD, aMatrix, 100, 0, 0, 0, [], [], 1, 0, W_TOL, D_LOWER_BOUND, [], [] );
[ D2, ~ ]= updateD_v8_ipopt( 'identity', 'L1', gY, WResStruct.W, WResStruct.W0, iD, DTemplate, BlkDS.indMap, 0, 1e-32, 1, 200, 5e-3, [] );

param.OUTER_IT_NUM=60;
[ expRec0 ] = dictionaryLearning_ADMM_v6( gY, initVar, DTemplate, [], 0, 0, 0, aMatrix, BlkDS, [], snapPath, param );
 LPAry(1) = LP_DL_Poiss( 'identity', aMatrix, gY, expRec0.W, expRec0.W0, gD, 0, 0, 0, 1, [], [], [] );
[ D0, ~ ]= updateD_v8_ipopt( 'identity', 'L1', gY, expRec0.W, expRec0.W0, gD, DTemplate, BlkDS.indMap, 0, 1e-32, 1, 200, 5e-3, [] );
 LPAry(2) = LP_DL_Poiss( 'identity', aMatrix, gY, expRec0.W, expRec0.W0, D0, 0, 0, 0, 1, [], [], [] );


[ ~, etheta, ~ ] = estimateHypParam( expRec0.W, expRec0.D, DTemplate, BlkDS, 0 );
[ maxTheta, minTheta ] = estimateL2ThetaMaxMin( gY, expRec0.D, expRec0.W, expRec0.W0, 'identity' );
[ expRec1 ] = dictionaryLearning_ADMM_v6( gY, initVar, DTemplate, [], 0, etheta*3.125, 0, aMatrix, BlkDS, [], snapPath, param );
[ D3, ~ ]= updateD_v8_ipopt( 'identity', 'L1', gY, expRec1.W, expRec1.W0, iD, DTemplate, BlkDS.indMap, 0, 1e-32, 1, 200, 5e-3, [] );

for i = 1:3
    figure;
    subplot(1,3, 1); imagesc(reshape(expRec1.W(i,:), 20, 20 ) ); colorbar; 
    subplot(1,3, 2); imagesc(reshape(gW(i,:), 20, 20 ) ); colorbar;
    subplot(1,3 ,3); imagesc(reshape(abs(expRec1.W(i,:)-gW(i,:)), 20, 20 ) ); colorbar;
end
for i = 1:3
    figure;
    subplot(1,3, 1); imagesc(reshape(expRec0.W(i,:), 20, 20 ) ); colorbar; 
    subplot(1,3, 2); imagesc(reshape(gW(i,:), 20, 20 ) ); colorbar;
    subplot(1,3 ,3); imagesc(reshape(abs(expRec0.W(i,:)-gW(i,:)), 20, 20 ) ); colorbar;
end
save('test_fusion_res.mat');

%% draw bar plot
[sg, idx] = sort( gD(:), 'descend');
outMat = zeros( length(gD(:)), 4);
outMat(:,1) = sg;
[ mfitD, Matching ] = HungarianArrange( expRec1.D, gD );
outMat(:,2) = mfitD(idx);
[ mfitD, Matching ] = HungarianArrange( expRec0.D, gD );
outMat(:,3) = mfitD(idx);
figure;bar(outMat)
legend('ground truth', 'DL with the fusion terms', 'DL w/o the fusion terms' );
xlabel('dictionary entries (sorted by intensity)', 'FontSize', 20, 'Fontweight', 'bold' );
ylabel('intensity', 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );

