gD = zeros(3,2); gD(1,1)=1; gD(2,1) = 1; gD(:,1) = gD(:, 1) / norm(gD(:,1), 2 ); gD(3,2) = 1;
DTemplate = zeros(3, 3); DTemplate(1:2, 1) = 1; DTemplate(1:2, 2) = 1; DTemplate(3, 3) = 1;
SLEN = 3; MLEN = 2; IHEIGHT = 20; IWIDTH = 20;
[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'diffusion', 1, [], 2, 1);
[ gY ] = genY_Poisson( gD, gW, gW0 );
BlkDS = conBLKDS( gY );

lambda = 1e-6; theta = 1e-6; phi = 1e-6;
aMatrix = ones( IHEIGHT, IWIDTH );
snapPath = 'temp.mat';
%setting up parameters
param = [];
param.OUTER_IT_NUM = 100; %number of  outer loop
param.ADMM_IT_NUM = 100; %number of ADMM loop
param.UP_D_IT_NUM = 200; %number of updating D loop
param.HES_FLAG = 0; %whether using real Hessian when updating D
param.CLUSTER_NUM = 12; %number of cluster uses, usually 12
param.SAVE_TEM_PERIOD = 1; %the perioud of iteration to save a temporary experiment result in snapPath
param.INNER_IT_RED_FLAG = 0; %whether reduce the ADMM_IT_NUM and UP_D_IT_NUM when in early outer loop
param.LATE_UPDATE_FLAG = 0; %whether using late update on dictionary elements
param.LATE_UPDATE_PERCENT = 0.2; %Under late update scenario, the percentage of dictionary elements used in the whole update process
param.LATE_UPDATE_INTTHRES = 0.8; %Under late update scenario, the percentage of intesity be explained used in the whole update process
param.CLUSTER_NAME = 'local'; %usually this is the default name
param.INIT_D = 'NNMF'; %the method of dictionary initialization
param.D_HIST_PATH = 'DHist_expRec.mat';
param.W_HIST_PATH = 'WHist_expRec.mat';

initVar = [];
initVar.outD = abs(randn(3,3) );
initVar.outD(DTemplate==0) = 0;
for i = 1:3
    initVar.outD = initVar.outD(:, i) / norm(initVar.outD(:, i) );
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


[ expRec ] = dictionaryLearning_ADMM_v4( gY, initVar, DTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( expRec.outW(i,:),20,20) ); colorbar;
end
[ elambda, ~, ephi ] = estimateHypParam( expRec.outW, expRec.outD, DTemplate, BlkDS );

[ expRec2 ] = dictionaryLearning_ADMM_v4( gY, initVar, DTemplate, [], elambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( expRec2.outW(i,:),20,20) ); colorbar;
end
[ ~, etheta, ~ ] = estimateHypParam( expRec2.outW, expRec2.outD, DTemplate, BlkDS );
[ expRec3 ] = dictionaryLearning_ADMM_v4( gY, initVar, DTemplate, [], elambda, etheta, dphi, aMatrix, BlkDS, [], snapPath, param );


%% both expRec4 and expRec5 work successfully deconvolute data
% these figures are generated from expRec6
load('DHist_expRec.mat');
for i = 1:length(DHistCell)
    if isempty( DHistCell{i} )
        break
    end
end
tIT_NUM = i-1;

outMat = zeros( tIT_NUM, 4 );
for i = 1:tIT_NUM
    outMat(i,1:2) = DHistCell{i}(1:2,1);
    outMat(i,3:4) = DHistCell{i}(1:2,2);
end
figure; bar(outMat)
legend('peak 1 in element 1', 'peak 2 in element 1', 'peak 1 in element 2', 'peak 2 in element 2')
title( 'Dictionary elements change over iterations' );
xlabel( '# iteration' );
ylabel( 'value' );
load('WHist_expRec.mat');
outMat2 = zeros( tIT_NUM, 2 );
for i = 1:tIT_NUM
    outMat2(i, 1) = max( WHistCell{i}(1,:) );
    outMat2(i, 2) = max( WHistCell{i}(2,:) );
end
figure; bar(outMat2);
legend('max weight of element 1', 'max weight of element 2', 'peak 1 in element 2', 'peak 2 in element 2')
title( 'maximun weights of elements change over iterations' );
xlabel( '# iteration' );
ylabel( 'value' );