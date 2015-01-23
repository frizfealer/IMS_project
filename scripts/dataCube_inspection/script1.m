%raster order
[sLen, hei, wid] = size(dataCube);
for i = 1:BlkDS.blkNum
    idx = BlkDS.B2GMap{i};
    curStart = 1;
    time = 1;
    rastIdx = zeros( length(idx), 1 );
    for j = 2:length(idx)
        if idx(j)-idx(j-1) == 1
            curEnd = j;
            if j == length(idx)
                if mod( time, 2 ) == 1
                    rastIdx(curStart:curEnd) = idx(curStart:curEnd);
                else
                    rastIdx(curStart:curEnd) = idx(curEnd:-1:curStart);
                end
            end
        elseif idx(j)-idx(j-1) > 1
            if mod( time, 2 ) == 1
                rastIdx(curStart:curEnd) = idx(curStart:curEnd);
            else
                rastIdx(curStart:curEnd) = idx(curEnd:-1:curStart);
            end
            curStart = j;
            time = time + 1;
        end
    end
    rastCell{i} = rastIdx;
end


totalIonMap = sum( dataCube, 1 );
totalIonMap = reshape( totalIonMap, hei, wid );
%i = 1
indMap = zeros( hei, wid);
% indMap( rastCell{i} ) = 1;
i=1;
for j = 1:8
    if j ==8
        sampleRange = (1+(j-1)*100):(length(rastCell{i}));
    else
        sampleRange = (1+(j-1)*100):(j*100);
    end
    indMap(rastCell{i}(sampleRange)) = totalIonMap(rastCell{i}(sampleRange));
    figure; subplot(2, 1, 1); imagesc(indMap); title( 'total count of each grid cells' );
    indMap(rastCell{i}(sampleRange)) = zeros(length(sampleRange), 1); colormap( bone ); colorbar;
    
    ins = dataCube(:, rastCell{i}(sampleRange));
    nnn = [];
    totalNum = 0;
    for k = 1:length(sampleRange)
        idx = find( ins(:, k) > 50 );
        nnn{k} = ins(idx, k);
        totalNum = totalNum + length(idx);
    end
    newIns = zeros(totalNum, 1);
    newStrVec = zeros(totalNum, 1);
    curNum = 1;
    for k = 1:length(sampleRange)
        curLen = length(nnn{k});
        newIns(curNum:(curNum+curLen-1)) = nnn{k};
        newStrVec(curNum:(curNum+curLen-1)) = repmat( k, curLen, 1);
        curNum = curNum + curLen;
    end
    %         newStrVec = num2str( newStrVec );
    [I,J] = ind2sub( [hei wid], rastCell{i}(sampleRange) );
    apop = ', ';
    leftS = '(';
    rightS = ')';
    I = num2str(I); J = num2str(J); labelVec = [repmat(leftS, length(I), 1) I repmat(apop, length(I), 1) J repmat(rightS, length(I), 1) ];
    subplot(2, 1, 2); boxplot( newIns, newStrVec, 'plotstyle', 'compact', 'label', labelVec ); title( 'boxplot of each m/z channels >= 50' );
     %         subplot(2, 1, 2); boxplot( ins, 'plotstyle', 'compact' );
end

i=2;
for j = 1:7
    if j ==7
        sampleRange = (1+(j-1)*100):(length(rastCell{i}));
    else
        sampleRange = (1+(j-1)*100):(j*100);
    end
    indMap(rastCell{i}(sampleRange)) = totalIonMap(rastCell{i}(sampleRange));
    figure; subplot(2, 1, 1); imagesc(indMap); title( 'total count of each grid cells' );
    indMap(rastCell{i}(sampleRange)) = zeros(length(sampleRange), 1); colormap( bone ); colorbar;
    
    ins = dataCube(:, rastCell{i}(sampleRange));
    nnn = [];
    totalNum = 0;
    for k = 1:length(sampleRange)
        idx = find( ins(:, k) > 50 );
        nnn{k} = ins(idx, k);
        totalNum = totalNum + length(idx);
    end
    newIns = zeros(totalNum, 1);
    newStrVec = zeros(totalNum, 1);
    curNum = 1;
    for k = 1:length(sampleRange)
        curLen = length(nnn{k});
        newIns(curNum:(curNum+curLen-1)) = nnn{k};
        newStrVec(curNum:(curNum+curLen-1)) = repmat( k, curLen, 1);
        curNum = curNum + curLen;
    end
    %         newStrVec = num2str( newStrVec );
    [I,J] = ind2sub( [hei wid], rastCell{i}(sampleRange) );
    apop = ', ';
    leftS = '(';
    rightS = ')';
    I = num2str(I); J = num2str(J); labelVec = [repmat(leftS, length(I), 1) I repmat(apop, length(I), 1) J repmat(rightS, length(I), 1) ];
    subplot(2, 1, 2); boxplot( newIns, newStrVec, 'plotstyle', 'compact', 'label', labelVec ); title( 'boxplot of each m/z channels >= 50' );
     %         subplot(2, 1, 2); boxplot( ins, 'plotstyle', 'compact' );
end

figure; subplot(2, 3, 1); plot( dataCube(:, 36, 42) ); title( 'spectrum at location (x=42, y=36)' );
subplot(2, 3, 2); plot( dataCube(:, 36, 43) ); title( 'spectrum at location (x=43, y=36)' );
subplot(2, 3, 3); plot( dataCube(:, 36, 44) ); title( 'spectrum at location (x=44, y=36)' );
subplot(2, 3, 4); plot( dataCube(:, 37, 42) ); title( 'spectrum at location (x=42, y=37)' );
subplot(2, 3, 5); plot( dataCube(:, 37, 43) ); title( 'spectrum at location (x=43, y=37)' );
subplot(2, 3, 6); plot( dataCube(:, 37, 44) ); title( 'spectrum at location (x=44, y=37)' );

figure; subplot(1,3,1); imagesc( reshape( sum(dataCube(313,:), 1), 38, 60 ) );colorbar; colormap(bone); title( ['signals distribution of m/z = ', num2str(mzAxis(313))] );
subplot(1,3,2); imagesc( reshape( sum(dataCube(314,:), 1), 38, 60 ) );colorbar; colormap(bone); title( ['signals distribution of m/z = ', num2str(mzAxis(314))] );
subplot(1,3,3); imagesc( reshape( sum(dataCube(313:314,:), 1), 38, 60 ) );colorbar; colormap(bone); title( ['signals distribution of m/z = ', num2str(mzAxis(313)) ' and m/z = ' num2str(mzAxis(314))] );

%% boxplot for each m/z values
for j = 1:4
    if j == 4
        mzRange = (1+(j-1)*100):( size(dataCube, 1) );
    else
        mzRange = (1+(j-1)*100):(j*100);
    end
    ins = dataCube(:, BlkDS.indMap == 1);
    ins = ins';
    tmp = mzAxis;
    nnn = [];
    totalNum = 0;
    for i = mzRange
        idx = find( ins(:, i) ~= 0 );
        nnn{i} = ins(idx, i);
        totalNum = totalNum + length(idx);
    end
    newIns = zeros(totalNum, 1);
    newStrVec = zeros(totalNum, 1);
    curNum = 1;
    for i =mzRange
        curLen = length(nnn{i});
        newIns(curNum:(curNum+curLen-1)) = nnn{i};
        newStrVec(curNum:(curNum+curLen-1), 1) = repmat( tmp(i), curLen, 1);
        curNum = curNum + curLen;
    end
    newStrVec = num2str( newStrVec ); newStrVec = newStrVec(:, 1:5);
    strVec = num2str(tmp); strVec = strVec(:, 1:5);
    figure; subplot(2, 1, 1); boxplot(newIns, newStrVec, 'plotstyle', 'compact'); title( 'boxplot of each m/z channels without zeros' );
    subplot(2, 1, 2); boxplot( ins(:, mzRange ), 'plotstyle', 'compact', 'label', strVec(mzRange, :) ); title( 'boxplot of each m/z channels with zeros' );
end

%% boxplot for different samples (locations)
for i = 1:8
    if i ~= 8
        figure; boxplot( ins(:, (1+(i-1)*50):(i*50) ), 'plotstyle', 'compact', 'label', strVec((1+(i-1)*50):(i*50),:) );
        ylim( [0, 9385] );
    else
        figure; boxplot( ins(:, (1+(i-1)*50):(size(ins, 2)) ), 'plotstyle', 'compact', 'label', strVec((1+(i-1)*50):(size(ins, 2)),:) );
        ylim( [0, 9385] );
    end
end