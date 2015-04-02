function [ aMatrix ] = leaveoutCBPatternsData_v2( BlkDS, percentage )
%--------------------------------------------------------------------------
% leaveoutCBPatternsData_v2: leave out Checker-Board pattern data version 2
%--------------------------------------------------------------------------
% DESCRIPTION:
%   create a indicator matrix of [height width] and with percentage of
%   leave-out data. The grids that are leaved out are 0; otherwise 1.
%
% INPUT ARGUMENTS:
%   BlkDS, block structure from conBLKDS
%   percentage, percentage of grids to be leaved out.
% OUTPUT ARGUMENTS:
%   aMatrix: size [h w], an indicator matrix
[hei, wid] = size(BlkDS.indMap);
aMatrix = zeros(hei, wid);
aMatrix(BlkDS.indMap==1) = 1;
nLenVec = zeros( BlkDS.blkNum, 1 );
for i = 1:BlkDS.blkNum
    nLenVec(i) = length(BlkDS.B2GMap{i});
    loLen = round( nLenVec(i)*percentage/100 );
    [I,J] = ind2sub( size(aMatrix), BlkDS.B2GMap{i} );
    heiRange = max(I) - min(I) + 1; widRange = max(J)-min(J) + 1;
    if heiRange < widRange
        ratio = round( widRange / heiRange );
        for intW= ratio*3:-1:3
            wNum = round( widRange / intW );
            hNum = round( loLen / wNum );
            for intH = 5:-1:3
                hRange = hNum * intH;
                if hRange < heiRange
                    break;
                end
            end
        end
    else
        ratio = round( heiRange / widRange );
        for intH= ratio*3:-1:3
            hNum = round( heiRange / intH );
            wNum = round( loLen / hNum );
            for intW = 5:-1:3
                wRange = wNum * intW;
                if wRange < widRange
                    break;
                end
            end
        end
    end
    loI = min(I):intH:max(I);
    loJ = min(J):intW:max(J);
    for k = loI
        for r = loJ
            aMatrix( k, r) = 0;
        end
    end
end

end

