function [ testRes, fig_Handle ] = data_sanity_check( inDCPath )
%inDCPath: datacube Path
[~,name,~] = fileparts( inDCPath );
load( inDCPath );
BlkDS = conBLKDS( dataCube );
%% check algorithm's correctness
insMap = zeros( size( dataCube, 2 ), size( dataCube, 3 ) );
for i = 1:size(posInfo, 1)
    insMap(posInfo(i,2),posInfo(i,1))=1;
end
% assert( isempty( find(insMap-BlkDS.indMap~=0, 1) ) );
fprintf( 'len of parameters = %d, len of data = %d\n', length(find(insMap==1)), length(find(BlkDS.indMap==1)) );
BLK_NUM = BlkDS.blkNum;

%% plot mean and total ion plot
[HEI, WID] = size( BlkDS.indMap );
ins = mean( dataCube, 1 ); ins = reshape( ins, HEI, WID );
fig_Handle = figure; subplot(2, 2, 1); imagesc(ins); colorbar; title( 'Mean Ion plot' );
ins = sum( dataCube, 1 ); ins = reshape( ins, HEI, WID );
subplot(2, 2, 2); imagesc(ins); colorbar; title( 'Total Ion plot' );
locInfo = cell( BLK_NUM, 1 );
for i = 1:BLK_NUM
    [locInfo{i}.I, locInfo{i}.J] = ind2sub( [HEI, WID], BlkDS.B2GMap{i} );
end

MarkerTable={'+','o','x','*', 's', 'd'};
ColorSpecTable=get(gca,'ColorOrder') ; 
usedColor = randperm( size(ColorSpecTable, 1) ); 
usedColor = usedColor(1:BLK_NUM); usedColor = ColorSpecTable(usedColor, :);

%% for the sample which should not be horizontally discriminant...
%test1: group testing
subplot(2, 2, 3);imagesc(ins); title( 'Two-sample t-test' );
for j = 1:BLK_NUM
    cI = locInfo{j}.I;
    cJ = locInfo{j}.J;
    uniCJ = unique( cJ );
    X = uniCJ';
    uY = zeros( 1, length(uniCJ) );
    lY = zeros( 1, length(uniCJ) );
    cnt = 1;
    for i = uniCJ'
        %fprintf('%d\n', i);
        cIdx = find(cJ==i);
        minI = min(cI(cIdx));
        maxI = max(cI(cIdx));
        uY(cnt) = minI;
        lY(cnt) = maxI;
        cnt = cnt + 1;
    end
    if j == 1
    [fillhandle,msg]=jbfill(X,uY,lY,usedColor(j, :),usedColor(j, :),1,1);
    else
    [fillhandle,msg]=jbfill(X,uY,lY,usedColor(j, :),usedColor(j, :),1,1);    
    end
end

h = zeros( BLK_NUM, 1 ); p = zeros( BLK_NUM, 1 );
textAry = cell( BLK_NUM, 1 );
for i = 1:BLK_NUM
    cI = locInfo{i}.I; cJ = locInfo{i}.J;
    boundary = floor( mean( min(cI):max(cI) ) );
    g1I = min(cI):boundary;
    g2I = (boundary+1):max(cI);
    g1Idx = find( ismember( cI, g1I ) == 1 );
    g2Idx = find( ismember( cI, g2I ) == 1 );
    g1Idx = sub2ind( [HEI, WID], cI(g1Idx), cJ(g1Idx) );
    g2Idx = sub2ind( [HEI, WID], cI(g2Idx), cJ(g2Idx) );
    [h(i),p(i)] = ttest2( ins(g1Idx), ins(g2Idx), 'Alpha', 0.01 );
    fprintf( 'h = %d, p-value = %g\n', h(i), p(i) );
    textAry{i} = sprintf( 'h = %d, p-value = %g\n',h(i), p(i) );
end
legend(textAry)
testRes{1} = [h p];
%% test2: line slope testing
subplot(2,2,4); title( 'Lint fitting test' );
slopeVec = zeros( BLK_NUM, 2 );
for i = 1:BLK_NUM
    cI = locInfo{i}.I; cJ = locInfo{i}.J;
    linVec = zeros( length(cI), 3 );
    for z = 1:length(cI)
        linVec(z, 1) = cI(z);
        linVec(z, 3) = ins(cI(z), cJ(z));
    end
    linVec(:, 2) = 1;
    slope = linVec(:,1:2)\linVec(:,3); slopeVec(i, :) = slope;
    cMarker = MarkerTable( mod(i, 6) + 1 );
    if i == 1
        scatter( linVec(:, 1), linVec(:, 3),  [], usedColor(i, :), cMarker{1} ); 
        hold on;
    else
        cH = plot( linVec(:, 1), linVec(:, 3), 'Color', usedColor(i, :), 'LineStyle', 'none' ); 
        set( cH,'Marker',cMarker{1});
    end
    res = linVec(:,3)-linVec(:,1:2)*slope;
    res = sum(res.^2)/length(res);
    fprintf( 'slope = %g, offset = %g res = %g\n',slope(1), slope(2), res );
    textAry{i} = sprintf( 'slope = %.3g, offset = %.3g res = %.3g\n',slope(1), slope(2), res );
end
lsline
testRes{2} = slopeVec;
hwplotprep;
print( fig_Handle, [name, '.pdf'], '-dpdf' );
close( fig_Handle );
% legend(textAry);