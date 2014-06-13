meanWAry = zeros(31,1);
for i = 1:31
    meanWAry(i) = mean(outW(i,:));
end
[val, idx] = sort(meanWAry, 'descend');

outW = full(outW);
outD = full(outD);
for j = 1:31
    figure;
    h1 = subplot(1,2,1);
    possIdx = find(mDTemplate(:,idx(j))~=0);
    mzRange = (mDMSIndex(possIdx(1),idx(j))-15):(mDMSIndex(possIdx(end),idx(j))+15);
    h=gca;
    plot(h, mzRange,zeros(length(mzRange), 1));
    for i = 1:length(possIdx)
        line( [mDMSIndex(possIdx(i),idx(j)) mDMSIndex(possIdx(i),idx(j))], [0 outD(possIdx(i), idx(j) )], 'LineWidth',2 );
        if outD(possIdx(i), idx(j) ) ~= 0
            if outD(possIdx(i), idx(j) ) >= 0.6
                annoText = ['[', deblank(mDIonName{idx(j)}(i,:)) ,']'];
                dist = length(annoText)/3.1*0.1;
                text(mDMSIndex(possIdx(i),idx(j))-3,max(outD(possIdx(i), idx(j) ))-dist,['[', deblank(mDIonName{idx(j)}(i,:)) ,']'], 'FontSize', 9, 'Rotation', 90);
            else
                text(mDMSIndex(possIdx(i),idx(j)),max(outD(possIdx(i), idx(j) ))+0.02,['[', deblank(mDIonName{idx(j)}(i,:)) ,']'], 'FontSize', 9, 'Rotation', 90);
            end
        end
    end
    title( ['Dictionary element: ', num2str(idx(j))]);
    xlabel( 'mzValue' ); ylabel( 'intensity' );
%     set(gca, 'ytick', sort(y(y~=0)));
    % set(gca, 'xtick', x(label(1):label(end)));
    h2 = subplot(1,2,2, 'align');
    temp1 = reshape(outW(idx(j),:,:),54, 88);
    imagesc(temp1);colorbar;
    title( ['The distribution in the space of Dictionary element: ', num2str(idx(j))]);
    xlabel( 'width' ); ylabel( 'height' );
    axis([h1 h2],'square')
    hwplotprep
    print ( '-dpdf', ['top_', num2str(j), '_DictionaryElement'] );
end

for i = 1:225
    possIdx = find(mDTemplate(:, i));
    MSVal = mDMSIndex(possIdx, i);
%     for j=1:length(MSVal)
%         fprintf('%d ', MSVal(j));
%     end
%     fprintf('\n');
    if(intersect([884, 898, 912, 926], MSVal))
        fprintf('%d\n', i);
    end
end

% generate name
W_name_cel = cell(14, 1);
for j = 1:10
    curIonType = mDIonName{idx(j)};
    possIdx = find(mDTemplate(:,idx(j))~=0);
    for i = 1:size(curIonType, 1)
        name = deblank( curIonType(i,:) );
        val = outD(possIdx(i), idx(j));
        if val >= 1e-6
            if strcmp(name, 'M+H') == 1
                W_name_cel{j} = [num2str( idx(j) ),': M+H=', num2str( mDMSIndex(possIdx(i),idx(j)) )];
                break;
            elseif strcmp(name, 'M+K') == 1
                W_name_cel{j} = [num2str( idx(j) ),': M+K=', num2str( mDMSIndex(possIdx(i),idx(j)) )];
                break;
            elseif strcmp(name, 'M+Na') == 1
                W_name_cel{j} = [num2str( idx(j) ),': M+Na=', num2str( mDMSIndex(possIdx(i),idx(j)) )];
                break;
            end
        end
    end
end
cnt = 11;
for j = [11 19 21 23]
    curIonType = mDIonName{idx(j)};
    possIdx = find(mDTemplate(:,idx(j))~=0);
    for i = 1:size(curIonType, 1)
        name = deblank( curIonType(i,:) );
        val = outD(possIdx(i), idx(j));
        if val >= 1e-6
            if strcmp(name, 'M+H') == 1
                W_name_cel{cnt} = [num2str( idx(j) ),': M+H=', num2str( mDMSIndex(possIdx(i),idx(j)) )];
                break;
            elseif strcmp(name, 'M+K') == 1
                W_name_cel{cnt} = [num2str( idx(j) ),': M+K=', num2str( mDMSIndex(possIdx(i),idx(j)) )];
                break;
            elseif strcmp(name, 'M+Na') == 1
                W_name_cel{cnt} = [num2str( idx(j) ),': M+Na=', num2str( mDMSIndex(possIdx(i),idx(j)) )];
                break;
            end
        end
    end
    cnt = cnt + 1;
end