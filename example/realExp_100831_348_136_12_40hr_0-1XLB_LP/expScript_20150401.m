% An experiment script for the experiment: 100831_348_136_12_40hr_0-1XLB_LP
%% set up file path
projectPath = 'D:\Users\YeuChern\GitHub\IMS_project';
inputFilePath = [ projectPath, '\', 'example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12,40hr_0-1XLB_LP.mzML_dc.mat' ];
basicVarPath = [ projectPath, '\', 'example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12_40hr_0_1XLB_LP_vars.mat' ];
%% prepare input variable for dictionary learning
exp_100831_348_136_12_40hr_0_1XLB_LP_basicSettings_v2( inputFilePath, basicVarPath );
%% Empirically test lambda, theta, phi's ranges.
%[ pRang ] = estimateHypParam( varFilePath, mode )
% lambda [0.2 20], theta [1e-4 1e-2], phi [1e5 1e8]
clear all;
projectPath = 'D:\Users\YeuChern\GitHub\IMS_project';
basicVarPath = [ projectPath, '\', 'example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12_40hr_0_1XLB_LP_vars.mat' ];
load( basicVarPath );
lambdaVec = logspace( log10(0.04), log10(40), 5 ); lambdaVec = [0 lambdaVec];
thetaVec = logspace( -4, -2, 5 ); thetaVec = [0 thetaVec];
phiVec = logspace( 3, 7, 5 ); phiVec = [0 phiVec];
save( basicVarPath ); clear all
%% run dictionary learning on all hyperparameters
%In current implementation, I use a function: 
%exp_100831_348_136_12_40hr_0_1XLB_LP_all_hypers_grid_search_v3
%And use a script:
%gen_hypers_grid_searchScript( 'D:\exp20150403_2', 'expScript_H_Pos', 'exp_100831_348_136_12_40hr_0_1XLB_LP_all_hypers_grid_search_v3', 216 )
%to generate scripts call the all_hypers_grid_search function.
%This way I can run all scripts in our clusters.
%% test result on leave out data
%choose the hyperparameters setting that has the best performace on the
%leave-out data
%Input:
%   resFolderPath: results folder that have all hyperparameters model.
%Output:
%    resMOLDLPath: has best MOLDL model and all hyperparameters testing
%    results
clear all;
resFolderPath = '/home/ycharn/all_hypers_grid_search_results_pos';
projectPath = '~/IMS_project-master';
resMOLDLPath = [ projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/all_hypers_grid_search_analysis_LP_20150403.mat' ];
basicVarPath = [ projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/100831_348_136_12_40hr_0_1XLB_LP_vars.mat' ];
[ lVec, tVec, pVec, valVec] = readResultsAndTest( resFolderPath, basicVarPath, 'identity', 1 );
save( resMOLDLPath, 'lVec', 'tVec', 'pVec', 'valVec' );

%% run sparse-NNMF on data
projectPath = '~/IMS_project-master';
basicVarPath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_vars.mat'];
outFilePath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LN/res_NNMF_20150403.mat'];
exp_100831_348_136_12_40hr_0_1XLB_LP_NNMF( basicVarPath, outFilePath );

%% run sparse-pLSA on data
projectPath = '~/IMS_project-master';
basicVarPath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_vars.mat'];
outFilePath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LN/res_pLSA_20150403.mat'];
exp_100831_348_136_12_40hr_0_1XLB_LP_pLSA( basicVarPath, outFilePath );

%% analyze the results
% MOLDL
resFolderPath = '/home/ycharn/all_hypers_grid_search_results_pos';
projectPath = '~/IMS_project-master';
resMOLDLPath = [ projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/all_hypers_grid_search_analysis_LP_20150403.mat' ];

load( resMOLDLPath );
tIdx = find(valVec==max(valVec));
listing = dir( resFolderPath );
disp(listing(tIdx+2).name); %exp_100831_348_136_12_40hr_0_1XLB_LP_res_v3_2_.mat
%load( [resFolderPath '/' listing(tIdx+2).name] );
%load the result directly
load( 'exp_100831_348_136_12_40hr_0_1XLB_LP_res_v3_2_.mat' );
W = expRec.W; D = expRec.D;
thresD = 1e-2; thresW = 1e-2;
[~, ~, molD, insW] = postProcessResults( D, W, thresD, thresW );
[~,ord] = sort( insW, 'descend');
%draw ten picture
% for i = 1:20
%     figure; imagesc(reshape(W(ord(i), :), 38, 59 ) ); colorbar;
% end

%compute surfactin reference spectrum
projectPath = '~/IMS_project-master';
SURFACTIN_INPUT_PATH = [projectPath '/' 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/30 ug LP.mzML_dc.mat'];
basicVarPath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/100831_348_136_12_40hr_0_1XLB_LP_vars.mat'];
MPPath = [projectPath, '/', 'example/molecule_profile_pos_v2.csv'];
load( SURFACTIN_INPUT_PATH );
load(basicVarPath);
figure; bar(dataVec, 0.01);

fileID = fopen(MPPath);
C = textscan(fileID, '%f %s', 'delimiter', {','});
fclose(fileID);
mpAry = C{1,1}; 
nameAry = C{1,2};

surfactinRef = [];
for i = 1:5
    surfactinRef(i).name = ['C' num2str(11+i)];
    surfactinRef(i).M = 994.3+(i-1)*14;
    surfactinRef(i).MIon = surfactinRef(i).M+mpAry;
    %take care the isotope that is not capture in ths referecne surfactin
    %data.
    if i == 4
        surfactinRef(i).MIon(1) = surfactinRef(i).MIon(1) - 1;
    end
    [mapIdx] = arrayfun( @(x)findingMathcedPeaks( mzAxis, x, 5e-4), surfactinRef(i).MIon );
    surfactinRef(i).mzIdx = mapIdx;
    surfactinRef(i).val = ones(length(mapIdx), 1)*-1;
    surfactinRef(i).val(mapIdx~=-1) = dataVec(mapIdx(mapIdx~=-1));
end
surfactinRef(4).MIon(1) = surfactinRef(4).MIon(1) + 1;

for i = 1:5
    [mapIdx] = arrayfun( @(x)findingMathcedPeaks( IMSD.mzAxis, x, 5e-4), surfactinRef(i).MIon);
    mapIMS = find(mapIdx~=-1);
    mapS = find( surfactinRef(i).mzIdx ~= -1 );
    iIdx = intersect(mapIMS, mapS);
    if ~isempty(iIdx)
        surfactinRef(i).existFlag = 1;
        refEle = zeros( size( IMSD.dataCube, 1), 1);
        %correct of peak matching.
        if i == 2
            mapIdx(3) = 422;
        end
        refEle(mapIdx(iIdx)) = surfactinRef(i).val(iIdx);
        refEle = refEle / norm( refEle );
        surfactinRef(i).refEle = refEle;
    else
        surfactinRef(i).existFlag = 0;
    end
end
%check cosine similarity of referece element to all
MOLDL_element = cell( 5, 1);
fitIdxVec = ones( 5, 1 )* -1;
for z = 1:5
    if surfactinRef(z).existFlag == 0
        continue;
    end
    refDele = surfactinRef(z).refEle;
    D = expRec.D;
    cosVec = zeros(size(D,2), 1);
    for i = 1:size(D,2)
        cosVec(i) = refDele'*D(:,i);
    end
    [~,fitIdxVec(z)]=max(cosVec);
%     if z == 4
%         fitIdxVec(z) = 624;
%     end
    MOLDL_element{z} = sparse( D(:, fitIdxVec(z)) );
end


for z = 1:5
    if surfactinRef(z).existFlag == 0
        continue;
    end
    cM = surfactinRef(z).MIon(1)-mpAry(1);
    CList = [];
    sList = [];
    for j = 1:length(SpeciesM2)
        iii = SpeciesM2{j};
        [mapIdx] = arrayfun( @(x)findingMathcedPeaks( iii, x, 5e-4), cM );
        if ~isempty( find(mapIdx~=-1, 1) )
            fprintf( 'z = %d, j = %d\n', z, j);
            [~,idx] = min(abs(cM-iii));
            CList = [CList, iii(idx)];
            sList = [sList, j];
        end
    end
    [~,idx ] = min(abs(cM-CList));
    fitIdxVec(z) = sList(idx);
    MOLDL_element{z} = sparse( D(:, fitIdxVec(z)) );
end

%analyze the result of spare-pLSA
projectPath = '~/IMS_project-master';
pLSAResPath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/res_pLSA_20150403.mat'];
load( pLSAResPath );
% [~, idx] = min(AICSPVec);
%choose the AIC that is the largest but with similar minima AIC
idx = 4;
D_pLSA = Pw_zSPVec{idx}; 
pLSA_element = cell( 5, 1);
for z = 1:5
    if surfactinRef(z).existFlag == 0
        continue;
    end
    refDele = surfactinRef(z).refEle;
    D = D_pLSA;
    cosVec = zeros(size(D,2), 1);
    for i = 1:size(D,2)
        cosVec(i) = refDele'*D(:,i);
    end
    [~,fitIdx]=max(cosVec);
    pLSA_element{z} = sparse( D(:, fitIdx) );
    pLSA_element{z}(pLSA_element{z}<max(pLSA_element{z})*1e-2) = 0;
    pLSA_element{z} = pLSA_element{z} / norm( pLSA_element{z} );
end

%analyze the result of sparse-NNMF
projectPath = '~/IMS_project-master';
NNMFResPath = [projectPath, '/', 'example/realExp_100831_348_136_12_40hr_0-1XLB_LP/res_NNMF_20150403.mat'];
load( NNMFResPath );
[~, idx] = min(ResSpVec);
NNMF_D = D_NNMF_SP_Cell{idx};
NNMF_element = cell( 5, 1 );
for z = 1:5
    if surfactinRef(z).existFlag == 0
        continue;
    end
    refDele = surfactinRef(z).refEle;
    cosVec = zeros(size(NNMF_D,2), 1);
    for i = 1:size(NNMF_D, 2)
        cosVec(i) = refDele'*NNMF_D(:,i);
    end
    [~,fitIdx]=max(cosVec);
    NNMF_element{z} = sparse( NNMF_D(:, fitIdx) );
     NNMF_element{z}(NNMF_element{z}<max(NNMF_element{z})*1e-2) = 0;
end


%plot all surfactin elements from different methods
for i = 2:4
figure; subplot(2, 2, 1);
plot(surfactinRef(i).refEle); 
subplot(2, 2, 2 );plot(MOLDL_element{i}); subplot(2, 2, 3 );plot(pLSA_element{i}); subplot(2, 2, 4 );plot(NNMF_element{i}); 
end

% plot the figures shown in the paper
outCell = cell( 3, 1 );
peakCell = cell( 3, 1 );
for i = 2:4
peakList = [];
idx = find( surfactinRef(i).refEle ~= 0 );
[~, vidx] = sort(surfactinRef(i).refEle(idx), 'descend' );
peakList = [peakList; idx(vidx)];
[~,ord]= sort( MOLDL_element{i}, 'descend' ); 
tmp = length(find(MOLDL_element{i} ~= 0));
if length(tmp) > 5
    tmp = 5;
end
peakList = [peakList; ord(1:tmp)];
[~,ord]= sort( pLSA_element{i}, 'descend' ); 
tmp = length(find(pLSA_element{i} ~= 0));
if tmp > 5
    tmp = 5;
end
peakList = [peakList; ord(1:tmp)];
[~,ord]= sort( NNMF_element{i}, 'descend' ); 
tmp = length(find(NNMF_element{i} ~= 0));
if tmp > 5
    tmp = 5;
end
peakList = [peakList; ord(1:tmp)];
peakList = unique(peakList);
peakCell{i} = peakList;
outMat = zeros(length(peakList),4);
for z = 1:length(peakList)
    outMat(z,:) = [surfactinRef(i).refEle(peakList(z)) MOLDL_element{i}(peakList(z)) pLSA_element{i}(peakList(z)) NNMF_element{i}(peakList(z)) ];
end
outCell{i} = outMat;
end

%plot
for i = 2:4
    outMat = outCell{i};
    peakList = peakCell{i};
    mzCVec = num2str(round(IMSD.mzAxis(peakList)*10)/10);
    figure;bar(outMat);
    set( gcf, 'Position', [1          25        1000        800] );
    legend( 'ground truth', 'MOLDL', 's-NNMF', 's-pLSA' );
    ylabel('intensity');
    ylabel( 'intesity', 'Fontsize', 20, 'Fontname','arial', 'fontweight', 'bold' );
    xlabel( 'm/z', 'Fontsize', 20, 'Fontname','arial', 'fontweight', 'bold' );
    set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
    set( gca, 'XTickLabel', mzCVec );
    xticklabel_rotate( [], 90, [], 'Fontsize', 20, 'Fontname','arial', 'fontweight','bold' );
    export_fig( ['surfactin-' surfactinRef(i).name, '.pdf'], '-transparent' );
end
