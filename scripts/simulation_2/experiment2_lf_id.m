gD = zeros(3,3); gD(1, 1) = 2; gD(2, 1) = 1; gD(2, 2) = 2; gD(3, 2) = 1; gD(1, 3) = 1; gD(3, 3) = 2;
for i = 1:3
gD = gD ./ norm(gD(:,i));
end
SLEN = 3; MLEN = 3; IHEIGHT = 20; IWIDTH = 20;
%[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'diffusion', 1, [], 2, 1);
%% construct W specific for experiment 2
gW = zeros( MLEN, IHEIGHT, IWIDTH );
gW(1, (IHEIGHT-4+1):IHEIGHT, 1:4) = abs(randn(4,4)*3e5);
midW = floor( IWIDTH / 2 );
gW(2, 1:4, (midW-2):(midW+1)) = abs(randn(4,4)*1e5);
gW(3, (IHEIGHT-2+1):IHEIGHT, (IWIDTH-2+1):IWIDTH ) = abs(randn(2,2)*3e6);
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

for i = 1:400
    for j =1:3
        if( gW(j, i) ~= 0)
            gW(j, i ) = (1000*abs(randn(1,1)))+gW(j,i);
        end
    end
end

figure;
for i = 1:3
subplot(1,3,i); imagesc(reshape(gW(i,:),20,20));colorbar;
end


[ gY ] = genY_Poisson( 'identity', gD, gW, gW0 );
BlkDS = conBLKDS( gY );

figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( gY(i, :), 20, 20 ) ); colorbar;
end

lambda = 1e-32; theta = 1e-32; phi = 1e-32;
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
save( 'experiment2_env_lf_id.mat' );


[ expRec0 ] = dictionaryLearning_ADMM_v6( gY, initVar, DTemplate, [], lambda, theta, phi, aMatrix, BlkDS, [], snapPath, param );
[ ~, etheta, ~ ] = estimateHypParam( expRec0.W, expRec0.D, DTemplate, BlkDS, 0 );
[ expRec1 ] = dictionaryLearning_ADMM_v6( gY, initVar, DTemplate, [], lambda, etheta*1e-2, phi, aMatrix, BlkDS, [], snapPath, param );

figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( expRec0.W(i,:),20,20) ); colorbar;
end
figure;
for i = 1:3
    subplot(1, 3, i); imagesc( reshape( expRec1.W(i,:),20,20) ); colorbar;
end
% [ elambda, ~, ephi ] = estimateHypParam( expRec0.outW, expRec0.outD, DTemplate, BlkDS );
% param.OUTER_IT_NUM = 60;
% [ expRec0_2 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], elambda, theta, ephi, aMatrix, BlkDS, [], snapPath, param );
% [ ~, etheta, ~ ] = estimateHypParam( expRec0_2.outW, expRec0_2.outD, DTemplate, BlkDS );
% 
% [ expRec1 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], elambda, etheta, ephi, aMatrix, BlkDS, [], snapPath, param );
% [ expRec3 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], lambda, etheta, phi, aMatrix, BlkDS, [], snapPath, param );
% [ ~, etheta, ~ ] = estimateHypParam( expRec0.outW, expRec0.outD, DTemplate, BlkDS );
% [ expRec4 ] = dictionaryLearning_ADMM_v5( gY, initVar, DTemplate, [], lambda, etheta, phi, aMatrix, BlkDS, [], snapPath, param );

%% compare with NNMF
opt = statset('MaxIter',100,'Display','final');
% [W, H] = nnmf( simData.gY(:, BlkDS.indMap==1), 20, 'algorithm', 'als', 'options', opt );
[W,H] = nnmf(gY(:, BlkDS.indMap==1), 3,'replicates',100,...
                   'options',opt,...
                   'algorithm','mult');
for i = 1:3
    W(:, i) = W(:, i) / norm( W(:, i) );
end
D_NNMF = W;


%% compare with pLSA
Par = []; Par.maxit = 500; Par.Leps = 1; Par.doplot = 0;
[Pw_z,Pd_z,Pz,Li] = pLSA_EM( gY(:, BlkDS.indMap==1), 3, Par );
    for j = 1:i
        Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
    end
D_pLSA = Pw_z;


%% draw bar plot
[sg, idx] = sort( gD(:), 'descend');
outMat = zeros( length(gD(:)), 4);
outMat(:,1) = sg;
[ mfitD, Matching ] = HungarianArrange( expRec0.D, gD );
outMat(:,2) = mfitD(idx);
[ mfitD, Matching ] = HungarianArrange( D_NNMF, gD );
outMat(:,3) = mfitD(idx);
[ mfitD, Matching ] = HungarianArrange( D_pLSA, gD );
outMat(:,4) = [mfitD(idx(1:9))];
figure;bar(outMat)
legend('ground truth', 'MOLDL', 'NNMF', 'pLSA');
xlabel('dictionary entries (sorted by intensity)', 'FontSize', 20, 'Fontweight', 'bold' );
ylabel('intensity', 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
save('simulation_res.mat');

