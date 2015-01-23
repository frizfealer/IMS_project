%% synthesize an example
gD = zeros(3,3); gD(1:2,1)=1; gD(2:3,2) = 1; gD([1],3) = 1;
for i = 1:3
    gD(:, i) = gD(:, i) / norm( gD(:, i) );
end
DTemplate = zeros(3, 3); DTemplate(1:2, 1) = 1; DTemplate(2:3, 2) = 1; DTemplate([1], 3) = 1;
SLEN = 3; MLEN = 3; IHEIGHT = 20; IWIDTH = 20;
[ gW, gW0, usedElement ] = synthesizeW( MLEN, IHEIGHT, IWIDTH, 'random', 1, 1, 10, 1);
[ gY ] = genY_Poisson( 'identity', gD, gW, gW0 );
BlkDS = conBLKDS( gY );
nullDTemplate = ones(3,3);
aMatrix = ones( IHEIGHT, IWIDTH );
%% test this example with dictionary learning without DTemplate
param = setDLParameters();
[ expRec1 ] = dictionaryLearning_ADMM_v6( gY, [], nullDTemplate, [], 1e-32, 1e-32, 1e-32, aMatrix, BlkDS, [], 'nullDTemp.mat', param );
%% test this example with dictionary learning with DTemplate
[ expRec2 ] = dictionaryLearning_ADMM_v6( gY, [], DTemplate, [], 1e-32, 1e-32, 1e-32, aMatrix, BlkDS, [], 'DTemp.mat', param );
%% draw bar plot
[sg, idx] = sort( gD(:), 'descend');
outMat = zeros( length(gD(:)), 4);
outMat(:,1) = sg;
[ mfitD, Matching ] = HungarianArrange( expRec2.D, gD );
outMat(:,2) = mfitD(idx);
[ mfitD, Matching ] = HungarianArrange( expRec1.D, gD );
outMat(:,3) = mfitD(idx);
figure;bar(outMat)
legend('ground truth', 'DL with the patterns', 'DL w/o the patterns' );
xlabel('dictionary entries (sorted by intensity)', 'FontSize', 20, 'Fontweight', 'bold' );
ylabel('intensity', 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
save('experiment_env.mat' );
