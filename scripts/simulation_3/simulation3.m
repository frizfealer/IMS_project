load('C:\Users\¦Ð±á\Google ¶³ºÝµwºÐ\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12_40hr_0_1XLB_LP_input.mat');
sDTemplate = pDTemplate; 
sDTemplate = sDTemplate(:,1:20);
tMZ = [];
for i = 1:size(sDTemplate, 1)
    if ~isempty( find( sDTemplate(i, :) ~= 0, 1 ) )
        tMZ = [tMZ i];
    end
end
sDTemplate = sDTemplate(tMZ, :);
save( 'simulation3_env.mat' );
[SLEN, MLEN] = size(sDTemplate);
HEIGHT = 30;
WIDTH = 30;
DOptions = [];
DOptions.coheMax = 0.5;
WOptions.type = 'diffusion';
WOptions.supPrec = 0.8;
[simData] = synthesizeData_Poisson( SLEN, MLEN, HEIGHT, WIDTH, sDTemplate, DOptions, WOptions, 1 );
BlkDS = conBLKDS(simData.gY);

%five-fold
% samNum = length(find(BlkDS.indMap(:)==1));
% samIdx = find(BlkDS.indMap(:)==1);
% foldNum = floor(samNum / 5);
% foldMap = zeros(30, 30);
% for i = 1:4
%     foldMap( samIdx( ((i-1)*foldNum+1):(foldNum*i)) ) = i;
% end
%five-fold only on 100 data point (two rectangle (y,x))
%(11, 5) (15, 14) and (16,21) (20,30)
%test = BlkDS.indMap==1;
%test(11:15, 5:14) = 0;
%test(16:20, 21:30) = 0;
%figure; imagesc(test);
%% using only partial data point of y (pY)
pY = zeros( size( simData.gY ) );
pY(:,11:15,5:14) = simData.gY(:, 11:15, 5:14 );
pY(:, 16:20, 21:30 ) = simData.gY(:, 16:20, 21:30 );
pBlkDS = conBLKDS(pY);
samNum = length(find(pBlkDS.indMap(:)==1));
samIdx = find(pBlkDS.indMap(:)==1);
foldNum = floor(samNum / 5);
foldMap = zeros(30, 30);
for i = 1:4
    foldMap( samIdx( ((i-1)*foldNum+1):(foldNum*i)) ) = i;
end
%figure;imagesc(foldMap);

% [ gridVec ] = genHypParamGrid( [9,9,9], 5 );
%% using special group of grids centering on elambda, etheta, ephi
[fVec, sVec, tVec] = meshgrid( 0:elambda:elambda*2, 0:etheta:etheta*2, 0:ephi:ephi*2 );
aMatrix = ones(30,30);
snapPath = 'simulation3_temp.mat';
param = setDLParameters();
param.LATE_UPDATE_FLAG = 0;
param.OUTER_IT_NUM = 30;
save( 'simulation3_env_1.mat' );

resVec = zeros(1:27, 2);
expTemp = [];
for i = 1:27
    lambda = fVec(i);
    theta = sVec(i);
    phi = tVec(i);
    vvv = zeros(5,1);
    for j = 1:5
        if j == 1
            param.OUTER_IT_NUM = 30;
            initVal = [];
        else
            param.OUTER_IT_NUM = 20;
            initVal = expRec;
        end
        aMatrix = ones(30, 30);
        aMatrix(foldMap(:)==j)=0;
        [ expRec ] = dictionaryLearning_ADMM_v5( pY, initVal, sDTemplate, [], lambda, theta, phi, aMatrix, pBlkDS, [], snapPath, param );
        expTemp{i,j} = expRec;
        save( 'resTemp.mat', 'expTemp', '-v7.3' );
        ttt=aMatrix==0;
        [ vvv(j) ] = validationOnTesting( ttt, pY, expRec.outW, expRec.outW0, expRec.outD );
    end
    resVec(i, 1) = mean(vvv);
    resVec(i, 2) = std(vvv);
    save( 'resVecTemp.mat', 'resVec' );
end

resVec_b = zeros(1:27, 2);
expTemp_b = [];
for i = 27:-1:1
    lambda = fVec(i);
    theta = sVec(i);
    phi = tVec(i);
    vvv = zeros(5,1);
    for j = 1:5
        if j == 1
            param.OUTER_IT_NUM = 30;
            initVal = [];
        else
            param.OUTER_IT_NUM = 20;
            initVal = expRec;
        end
        aMatrix = ones(30, 30);
        aMatrix(foldMap(:)==j)=0;
        [ expRec ] = dictionaryLearning_ADMM_v5( pY, initVal, sDTemplate, [], lambda, theta, phi, aMatrix, pBlkDS, [], snapPath, param );
        expTemp_b{i,j} = expRec;
        save( 'resTemp_b.mat', 'expTemp_b', '-v7.3' );
        ttt=aMatrix==0;
        [ vvv(j) ] = validationOnTesting( ttt, pY, expRec.outW, expRec.outW0, expRec.outD );
    end
    resVec_b(i, 1) = mean(vvv);
    resVec_b(i, 2) = std(vvv);
    save( 'resVecTemp_b.mat', 'resVec_b' );
end

