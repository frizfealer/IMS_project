BlkDS = conBLKDS( tDataCube );
BLK_NUM = BlkDS.blkNum;
[HEI, WID] = size( BlkDS.indMap );
ins = mean( tDataCube, 1 ); ins = reshape( ins, HEI, WID );
figure; subplot(1, 2, 1); imagesc(ins); colorbar; title( 'Mean Ion plot' );
ins = sum( tDataCube, 1 ); ins = reshape( ins, HEI, WID );
subplot(1, 2, 2); imagesc(ins); colorbar; title( 'Total Ion plot' );
locInfo = cell( BLK_NUM, 1 );
for i = 1:BLK_NUM
    [locInfo{i}.I, locInfo{i}.J] = ind2sub( [HEI, WID], BlkDS.B2GMap{i} );
end

%% for the sample which should not be horizontally discriminant...
%test1: group testing
h = zeros( BLK_NUM, 1 ); p = zeros( BLK_NUM, 1 );
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
    fprintf( 'grp %d: h = %d, p-value = %g\n', i, h(i), p(i) );
end
%test2: line slope testing
for i = 1:BLK_NUM
    cI = locInfo{i}.I; cJ = locInfo{i}.J;
    linVec = zeros( length(cI), 3 );
    for z = 1:length(cI)
        linVec(z, 1) = cI(z);
        linVec(z, 3) = ins(cI(z), cJ(z));
    end
    linVec(:, 2) = 1;
    slope = linVec(:,1:2)\linVec(:,3);
    res = linVec(:,3)-linVec(:,1:2)*slope;
    res = sum(res.^2)/length(res);
    fprintf( 'grp %d: slope = %g, offset = %g res = %g\n', i, slope(1), slope(2), res );
end
    