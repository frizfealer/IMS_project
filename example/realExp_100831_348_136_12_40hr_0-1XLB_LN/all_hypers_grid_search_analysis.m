function [] = all_hypers_grid_search_analysis( testingFlag )
%% analyze all hyperparameters grid search results
% INPUT_PATH = '~/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat';
% RESULT_FOLDER = '~/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/all_hypers_grid_search_results';
INPUT_PATH = 'D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat';
RESULT_FOLDER = 'D:\IMS_DATA\100831_348_136_12,40hr_0-1XLB_LN\all_hypers_grid_search_results';
%% test objetive function on all hyperparameters models
if testingFlag == 1
    [ lVec, tVec, pVec, valVec ] = readResultsAndTest( RESULT_FOLDER, INPUT_PATH, 'identity', 1 );
    save( 'all_hypers_grid_search_analysis_LN.mat', 'lVec', 'tVec', 'pVec', 'valVec' );
end
%% reload the results on any computers
load( 'all_hypers_grid_search_analysis_LN.mat' );
tIdx = find(valVec==max(valVec));
listing = dir( RESULT_FOLDER );
listing(tIdx+2).name
load( [RESULT_FOLDER '\' listing(tIdx+2).name] );
%tIdx =26
W = expRec.W; D = expRec.D;
thresD = 1e-2; thresW = 1e-2;
[~, ~, molD, insW] = postProcessResults( D, W, thresD, thresW );
[~,ord] = sort( insW, 'descend');
for i = 1:20
    figure; imagesc(reshape(W(ord(i), :), 38, 60 ) ); colorbar;
end

%% compute surfactin reference spectrum
SURFACTIN_INPUT_PATH = 'surfactin_RN.mzML_dc.mat';
load( SURFACTIN_INPUT_PATH );
load(INPUT_PATH);
rDataVec_RN(rDataVec_RN<max(rDataVec_RN)*1e-2) = 0;
figure; plot(rDataVec_RN);
[~,ord] = sort( rDataVec_RN, 'descend');
rMZIdx = sort(ord(1:5));
refDelVal = sqrt( rDataVec_RN( rMZIdx ).^2/sum(rDataVec_RN( rMZIdx).^2) );
refDele = zeros( size(dataCube, 1), 1 );
idx = find(refDelVal >= max(refDelVal)*1e-2);
refDelVal = refDelVal(idx);
rMZIdx = rMZIdx(idx);
%inspect the related peaks index in experiment's mzAxis
INS_THRES = 1e-2;
mzTarget = rMZAxis_RN(rMZIdx);
[ tIdxVec ] = check_related_peaks_index( mzTarget, mzAxis, W );
ins = zeros( size(dataCube, 1), 1 );
for i = 1:size(dataCube, 1)
    ins(i) = mean( dataCube(i, BlkDS.indMap == 1 ) );
end
maxIns_in_data = max(ins);
ins = zeros( length(tIdxVec), 1 );
for i = 1:length(tIdxVec)
    figure; imagesc( reshape(dataCube(tIdxVec(i),:),38,60));
    ins(i) = mean( dataCube(tIdxVec(i), BlkDS.indMap == 1 ) );
end
idx = find(ins>maxIns_in_data*INS_THRES);
tIdxVec = tIdxVec(idx);
refDele(tIdxVec) = refDelVal(idx); refDele = sparse(refDele);
%% check cosine similarity of referece element to all
D = expRec.D;
cosVec = zeros(size(D,2), 1);
for i = 1:size(D,2)
    cosVec(i) = refDele'*D(:,i);
end
[~,fitIdx]=max(cosVec);
MOLDL_element = sparse( D(:, fitIdx) );

%% load pLSA results
D_THRES = 1e-2;
load('D_pLSA_LL_AIC_LN.mat')
[~,best_pLSA] = min(AICVec(1:length(D_pLSA_Cell)));
pLSA_D = D_pLSA_Cell{best_pLSA};
cosVec = zeros(size(pLSA_D,2), 1);
for i = 1:size(pLSA_D,2)
    cosVec(i) = refDele'*pLSA_D(:,i);
end
[~,fitIdx]=max(cosVec);
pLSA_element = sparse(pLSA_D(:, fitIdx));
pLSA_element(pLSA_element<max(pLSA_element)*D_THRES) = 0;

%% load NNMF results
load('NNMF_LN.mat')
resVec = zeros( length(D_NNMF_Cell), 1 );
%from 8 to 150
for i = 8:length(resVec)
    resVec(i) = Res_Cell{i};
end
[~,best_NNMF] = min(resVec(8:end)); best_NNMF = best_NNMF+8;
NNMF_D = D_NNMF_Cell{best_NNMF};
cosVec = zeros(size(NNMF_D,2), 1);
for i = 1:size(NNMF_D,2)
    cosVec(i) = refDele'*NNMF_D(:,i);
end
[~,fitIdx]=max(cosVec);
NNMF_element = sparse(NNMF_D(:, fitIdx));
NNMF_element(NNMF_element<max(NNMF_element)*D_THRES) = 0;

%plot all surfactin elements from different methods
figure; subplot(2, 2, 1);
plot(refDele); 
subplot(2, 2, 2 );plot(MOLDL_element); subplot(2, 2, 3 );plot(pLSA_element); subplot(2, 2, 4 );plot(NNMF_element); 

%% plot the figures shown in the paper
peakList = [];
[~,ord]= sort( refDele, 'descend' ); 
peakList = [peakList; ord(1:length(tIdxVec))];
[~,ord]= sort( MOLDL_element, 'descend' ); 
tmp = length(find(MOLDL_element ~= 0));
if length(tmp) > 5
    tmp = 5;
end
peakList = [peakList; ord(1:tmp)];
[~,ord]= sort( pLSA_element, 'descend' ); 
tmp = length(find(pLSA_element ~= 0));
if length(tmp) > 5
    tmp = 5;
end
peakList = [peakList; ord(1:tmp)];
[~,ord]= sort( NNMF_element, 'descend' ); 
tmp = length(find(NNMF_element ~= 0));
if length(tmp) > 5
    tmp = 5;
end
peakList = [peakList; ord(1:tmp)];
peakList = unique(peakList);

outMat = zeros(length(peakList),4);
for i = 1:length(peakList)
    outMat(i,:) = [refDele(peakList(i)) MOLDL_element(peakList(i)) pLSA_element(peakList(i)) NNMF_element(peakList(i)) ];
end

%% plot
mzCVec = num2str(round(mzAxis(peakList)));
figure;bar(outMat);
legend( 'ground truth', 'DL', 'NNMF', 'pLSA' );
ylabel('intensity');
ylabel( 'intesity', 'Fontsize', 20, 'Fontname','arial' );
set( gca, 'FontSize', 20, 'Fontweight', 'bold' );
set( gca, 'XTickLabel', mzCVec );
xticklabel_rotate( [], 90, [], 'Fontsize', 20, 'Fontname','arial', 'fontweight','bold' );

end

