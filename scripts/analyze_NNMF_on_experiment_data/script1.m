%load input variables for negative mode data.
load( 'D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );

[~, pIdx1] = min(abs(mzAxis-1007.7));
[~, pIdx2] = min(abs(mzAxis-1021.8));
[~, pIdx3] = min(abs(mzAxis-1035.8));
% A possible error peak
[~, pIdx4] = min(abs(mzAxis-1057.6));
%% shown in the data
figure; subplot( 2, 2, 1); imagesc( reshape( dataCube(pIdx1, :), 38, 60 ) ); colorbar ;title('m/z = 1007.7');
subplot( 2, 2, 2 ); imagesc( reshape( dataCube(pIdx2, :), 38, 60 ) ); colorbar ;title('m/z = 1021.8');
subplot( 2, 2, 3 ); imagesc( reshape( dataCube(pIdx3, :), 38, 60 ) ); colorbar ;title('m/z = 1035.8');
subplot( 2, 2, 4 ); imagesc( reshape( dataCube(pIdx4, :), 38, 60 ) ); colorbar ;title('m/z = 1057.6');
export_fig four_peaks_in_data.pdf -transparent
%% get target compounds in dictionary 
targetM = [];
for i = 1:size( nDTemplate, 2)
    if nDTemplate(pIdx1, i) == 1
        targetM = [targetM, i];
    elseif nDTemplate(pIdx2, i) == 1
        targetM = [targetM, i];
    elseif nDTemplate(pIdx3, i) == 1
        targetM = [targetM, i];
    elseif nDTemplate(pIdx4, i) == 1
        targetM = [targetM, i];
    end
end
initD = nDTemplate(:, targetM);

targetM2 = [];
for i = 1:size( nDTemplate, 2)
    if nDTemplate(pIdx1, i) == 1 && nDTemplate(pIdx2, i) == 1 && nDTemplate(pIdx3, i) == 1
        targetM2 = [targetM2, i];
    end
end
initD2 = nDTemplate(:, targetM2);


%% test on NNMF
inY = dataCube(:, :);
D_NNMF_Cell = {};
Res_Cell = {};
for i = 8:size(inY, 1)
mNum = i;
rNum = 1e3;
[ D1, R1] = NNMF_als_wrapper( inY, [], BlkDS, mNum, rNum, true, 'default' );
D_NNMF_Cell{i} = sparse(D1); Res_Cell{i} = R1;
%run on lab machine
save( 'res_12152014.mat', 'D_NNMF_Cell', 'Res_Cell' );
% figure; subplot( 2, 2, 1); plot( D_NNMF_Cell{i}(pIdx1, :) ); title('m/z = 1007.7');
% subplot( 2, 2, 2 ); plot( D_NNMF_Cell{i}(pIdx2, :) ); title('m/z = 1021.8');
% subplot( 2, 2, 3 ); plot( D_NNMF_Cell{i}(pIdx3, :) ); title('m/z = 1035.8');
% subplot( 2, 2, 4 ); plot( D_NNMF_Cell{i}(pIdx4, :) ); title('m/z = 1057.6');
end

%% generage figure: residuals versus k-number
ResVec = zeros( length(Res_Cell), 1 );
for i = 1:length(Res_Cell)
    ResVec(i) = Res_Cell{i};
end
figure; plot(ResVec); xlabel( 'compound number' ); ylabel( 'residuals' );

%% evaluate the dictionary recovery based on three peaks
statVec = zeros( length(Res_Cell), 1 );
for i = 1:length(Res_Cell)
    [ resVec ] = evaluateRecovery( D_NNMF_Cell{i}, [pIdx1, pIdx2, pIdx3], setdiff(1:396, [pIdx1, pIdx2, pIdx3]) );
    statVec(i) = max(resVec);
end
[ ~, maxK ] = max(statVec);
[ resVec ] = evaluateRecovery( D_NNMF_Cell{maxK}, [pIdx1, pIdx2, pIdx3], setdiff(1:396, [pIdx1, pIdx2, pIdx3]) );
[ ~, maxE ] = max( resVec );

tM = maxE;
figure; plot( D_NNMF_Cell{maxK}(:, tM) ); xlabel( 'm/z value' ); ylabel( 'intensity' );
title( ['spectrum of element ', num2str(tM), ' in NNMF, compound number = ', num2str(maxK) ] );
xTickCell = {};
cnt = 1;
curXTick =  get(gca, 'XTickLabel'); curXTick = str2num(curXTick);
for i = 1:length(curXTick)
    if i == 1 || i == length(curXTick)
        xTickCell{i} = '';
    else
        xTickCell{i} = mzAxis((curXTick(i)));
    end
end
set(gca,'XTickLabel',xTickCell)
export_fig dict_element_captured_peaks.pdf -transparent
