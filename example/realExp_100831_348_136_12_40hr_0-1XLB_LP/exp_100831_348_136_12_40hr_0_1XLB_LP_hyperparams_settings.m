function [] = exp_100831_348_136_12_40hr_0_1XLB_LP_hyperparams_settings
%PROJECT_FOLDER_PATH: you should change to your own folder of the codes.
PROJECT_FOLDER_PATH = 'D:\Users\YeuChern\GitHub\IMS_project';
load([PROJECT_FOLDER_PATH '\' ...
    'example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12_40hr_0_1XLB_LP_vars.mat'] );

%% construct a subset of samples (InDa_set1, InDa_set2, sWA_set1, sWA_set2)
%find the cell grid with the largest signal.
ins = max( InDa.dataCube, [], 1 ); 
ins = reshape( ins, size(ins, 2), size(ins, 3) );
figure; imagesc(ins);
[i, j] = find( ins == max(ins(:)) );
tIdx = sub2ind([size(ins, 1) size(ins,2)], i, j);
%find the colony of this grid cell.
cBlk = InDa.BlkDS.G2BMap{tIdx};
%get the bound of the colony.
[bI, bJ] = ind2sub( [size(ins, 1) size(ins,2)], InDa.BlkDS.B2GMap{cBlk} );
%height equals to the whole height of this colony,
%width equals to 21.
cH = j + 10; cL = j - 10;
if cH > max(bJ)
    cDiff = cH - max(bJ) + 1;
    cH = max(bJ);
    cL = cL - cDiff;
    assert( cL >= min(bJ) );
elseif cL < min(bJ)
    cDiff = min(bJ) - cL + 1;
    cL = min(bJ);
    cH = cH + cDiff;
    assert( cH <= max(bJ) );
end
%construct the InDa subset 1 and the respective WA data structure.
hidRang = min(bI):max(bI);
widRang = cL:cH;
InDa_set1 = initInputData( InDa.dataCube(:, hidRang, widRang),...
    InDa.mzAxis, 0);
sWA_set1 = initWADMM( InDa_set1, size(sDI.DTemplate, 2) );
% 
usedCell = zeros(length(hidRang)*length(widRang), 1);
cnt = 1;
for i = 1:length(hidRang)
    for j = 1:length(widRang)
        usedCell(cnt) = find((bI == hidRang(i)) & (bJ == widRang(j)));
        cnt = cnt + 1;
    end
end
rBI = bI; rBI(usedCell) = []; rBJ = bJ; rBJ(usedCell) = [];
rBI = rBI(1:size(InDa_set1.dataCube, 2)*size(InDa_set1.dataCube, 3)-1);
rBJ = rBJ(1:size(InDa_set1.dataCube, 2)*size(InDa_set1.dataCube, 3)-1);
tmp = sub2ind( [size(ins, 1) size(ins,2)], rBI, rBJ );
tmp = InDa.dataCube(:, tmp); 
InDa_set2 = initInputData( tmp,...
    InDa.mzAxis, 0);
sWA_set2 = initWADMM( InDa_set2, size(sDI.DTemplate, 2) );

%% testing lambda's range
lambdaVec = linspace(  1, 5, 100 ); lambdaVec = [0 lambdaVec];
AStruc = cell( length(lambdaVec), 1 );
parfor i = 1:length(lambdaVec)
    fprintf( '%d\n', i );
    cWA = sWA_set1; cWA.lambda = lambdaVec(i); cWA.itNum = 100;
    tmp = updateW_ADMM_v7( InDa_set1, sDI.DTemplate, cWA, 0 );
    AStruc{i} = tmp.W;
end
for i = 1:10:101
    insW = AStruc{i}; insW = max( insW, [], 1 );
    figure; imagesc(reshape(insW, size(InDa_set1.dataCube, 2), size(InDa_set1.dataCube, 3)) ); 
    title(['l = ', num2str(lambdaVec(i))]);
end
%lambda should be 0, 1~3
lambdaVec = linspace( 1, 4, 7 ); lambdaVec = [0 lambdaVec];

%% testing theta's range
thetaVec = logspace( -5, -2, 100 ); thetaVec = [0 thetaVec];
ASruct = cell( length(thetaVec), 1 );
parfor i = 1:length(thetaVec)
    fprintf( '%d\n', i );
    cWA = sWA_set1; cWA.theta = thetaVec(i); cWA.itNum = 100;
    tmp = updateW_ADMM_v7( InDa_set1, sDI.DTemplate, cWA, 0 );    
    AStruc{i} = tmp.W;
end
for i = 1:10:101
    insW = AStruc{i}; insW = max( insW, [], 1 );
    figure; imagesc(reshape(insW, size(InDa_set1.dataCube, 2), size(InDa_set1.dataCube, 3) ) ); 
    title(['t = ', num2str(thetaVec(i))]);
end
%theta should be 0, 1.0000e-05~1.2328e-05
thetaVec = logspace( log10(1e-6), log10(1e-4), 7 ); thetaVec = [0 thetaVec];

%% testing phi's range
phiVec = logspace( 4, 8, 100 ); phiVec = [0 phiVec];
DStruc = cell( length(phiVec), 1 );
cWA = sWA_set1; 
tmpWA = updateW_ADMM_v7( InDa_set1, sDI.DTemplate, cWA, 0 );
parfor i = 1:length(phiVec)
    fprintf( '%d\n', i );
    cDI = sDI; cDI.phi = phiVec(i);
    [tmp, ~] = updateD_v9_ipopt( InDa_set1, ...,
        tmpWA.W, tmpWA.W0, cDI, 0 );
    DStruc{i} = tmp.D;
end
for i = 1:10:91
    tmp = DStruc{i}; tmp(tmp<=0.01) = 0;
    tmp2 = DStruc{i+10}; tmp2(tmp2<=0.01)=0;
    lenVec = zeros(length(sDI.DIonName), 1);
    normVec = zeros(length(sDI.DIonName), 1);
    for j = 1:length(sDI.DIonName)
        lenVec(j) = length(find(tmp(:, j) ~=0));
        normVec(j) = norm(tmp(:, j));
    end
    %figure; plot(lenVec);
    figure;imagesc(tmp-tmp2); colorbar;
    %figure;plot(normVec);
    title(['p  ', num2str(phiVec(i)), ' - p ', num2str(phiVec(i+10))]);
end
%we can check the difference between D to get the upper bound
%phi should be 0, 1e4~4e6
phiVec = logspace( log10(1e3), log10(4e6), 7 ); phiVec = [0 phiVec];
%% save variabe 
%InputData structure "InDa", "sDI", the start DI, "sWA", the start WA.
%"InDa_set1", "InDa_set2", "sWA_set1", "sWA_set2" for grid search and testing
%
OUTPUT_FILE_PATH = 'example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12_40hr_0_1XLB_LP_vars.mat';
save( [PROJECT_FOLDER_PATH '\' OUTPUT_FILE_PATH],  'InDa', 'sDI', 'sWA', ...
    'InDa_set1', 'sWA_set1', ...
    'InDa_set2', 'sWA_set2', ...
    'lambdaVec', 'thetaVec', 'phiVec' );

end