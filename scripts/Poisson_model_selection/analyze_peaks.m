function [ resCell ] = analyze_peaks( DTemplate, dataCube, mzAxis, tElement, outFileFlag )
%analyze_peaks Summary of this function goes here
%tElement: target dictionary elements, can be scalar or vector
%outFilePath, the path to output image file
ROW_PIC_NUM = 4;
IGNORE_INT = 1e-2;
[~, hei, wid] = size( dataCube );

for i = 1:length(tElement)
    idx = find( DTemplate(:,tElement(i))>IGNORE_INT );
    h = figure;
    eNum = length(idx);
    [rNum, cNum] = determineRC( eNum, ROW_PIC_NUM );
    for j = 1:eNum
        subplot( rNum, cNum, j );
        imagesc( reshape( dataCube( idx(j), :), hei, wid ) ); colorbar;
        title( ['m/z = ', num2str(mzAxis(idx(j)))] );
    end
    curCC = zeros( eNum, eNum );
    for j = 1:(eNum-1)
        for k = (j+1):eNum
            temp = corrcoef( dataCube(idx(j), :), dataCube(idx(k), :) );
            curCC(j, k) = temp(1, 2);
        end
    end
    if outFileFlag == 1
        hwplotprep;
        print( '-dpdf', [num2str(tElement(i)), '_different_peaks.pdf'] );
        close(h);
    end
    resCell{i} = curCC;
end

end

function [rowN, colN] = determineRC( tNum, row_pic_num )
    rowN = floor( (tNum-1) / row_pic_num ) + 1;
    colN = row_pic_num;
end

