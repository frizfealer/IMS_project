function [ mSpec, mMZAxis, binNumVec ] = binningSpectrum( spec, mzAxis )
    sLen = length(spec);
    ins = spec;
    [~,dOrd] = sort(ins, 'descend');
    tol = 0.5;
    dMZ = mzAxis(dOrd);
    mMZAxis = zeros( length(dMZ), 1 );
    mIndVec = zeros( length(dMZ), 1 );
    curPeakID = 1;
    binNumVec = [];
    while ~isempty(dMZ)
%         fprintf( '%d\n', curPeakID );
%         idx = ismember( mzAxis, dMZ ) ;
%         tmp = ins(idx);
%         if max(tmp) ~= ins( mzAxis == dMZ(1) )
%             keyboard();
%         end
        lMZ = dMZ(1);
        rMZ = dMZ( dMZ >= lMZ - tol );
        rMZ = rMZ( rMZ <= lMZ + tol );
        mIndVec( ismember( mzAxis, rMZ )==1 ) = curPeakID;
        mMZAxis(curPeakID) = lMZ;
        binNumVec(curPeakID) = length(rMZ);
        curPeakID = curPeakID + 1;
        rDOrd = find(ismember( dMZ, rMZ ));
%         fprintf( '%d, %d, %d\n', length(dMZ), length(rMZ), length(rDOrd) );
%         if length(rMZ) ~= length(rDOrd)
%             keyboard();
%         end
        dMZ(rDOrd) = [];
    end

    peakNum = max(mIndVec);
    tmp = mMZAxis;
    mMZAxis = tmp(1:peakNum);
    mMZAxis = sort( mMZAxis );
    mSpec = zeros( peakNum, 1);
    startPos = 1;
    curPeakID = 1;
    for j = 1:length(mIndVec)-1
        if mIndVec(j+1) ~= mIndVec(j)
            endPos = j;
            mSpec(curPeakID) = max(spec(startPos:endPos));
            startPos = j+1;
            curPeakID = curPeakID + 1;
        end
    end
    endPos = j + 1;
    mSpec(curPeakID) = max(spec(startPos:endPos));
end