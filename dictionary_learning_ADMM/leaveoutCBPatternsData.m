function [ aMatrix ] = leaveoutCBPatternsData( BlkDS, percentage )
%leaveoutCBPatternsData leave out Checker-Board pattern data
[hei, wid] = size(BlkDS.indMap);
aMatrix = zeros(hei, wid);
aMatrix(BlkDS.indMap==1) = 1;
nLenVec = zeros( BlkDS.blkNum, 1 );
for i = 1:BlkDS.blkNum
    nLenVec(i) = length(BlkDS.B2GMap{i});
    loLen = round( nLenVec(i)*percentage/100 );
    [I,J] = ind2sub( size(aMatrix), BlkDS.B2GMap{i} );
    heiRange = max(I) - min(I); widRange = max(J)-min(J);
    if heiRange < widRange
        startHei = min(I) + 2;
        endHei = max(I) - 2;
        hLen = endHei - startHei + 1;
        wLen = round( 2*loLen / hLen );
        startWid = ceil( (max(J)-min(J))/2 )+min(J) - floor(wLen/2);
        endWid = ceil( (max(J)-min(J))/2 )+min(J) + floor(wLen/2);
        assert( startWid >= min(J) ); assert( endWid <= max(J) );
    else
        startWid = min(J) + 2;
        endWid = max(J) - 2;
        wLen = endWid - startWid + 1;
        hLen = round( 2*loLen / wLen );
        startHei = ceil( (max(I)-min(I))/2 )+min(I) - floor(hLen/2);
        endHei = ceil( (max(I)-min(I))/2 )+min(I) + floor(hLen/2);
        assert( startHei >= min(I) ); assert( endHei <= max(I) );
    end
    cMat = checkerboard_1( hLen,  wLen );
    [I, J] = find( cMat == 1 );
    I = I + startHei;
    J = J + startWid;
    aMatrix( sub2ind( size(aMatrix), I, J) ) = 0;
end

end

