function [ mDataCube, mMZAxis, iii ] = binningDataCube( dataCube, mzAxis, BlkDS )
    [~, hei, wid] = size(dataCube);
    ins = dataCube(:, BlkDS.indMap == 1 );
    ins = sum( ins, 2 );
    [~,dOrd] = sort(ins, 'descend');
    tol = 0.5;
    dMZ = mzAxis(dOrd);
    cMZ = mzAxis;
    mMZAxis = zeros( length(dMZ), 1 );
    mIndVec = zeros( length(dMZ), 1 );
    curPeakID = 1;
    iii = [];
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
        iii(curPeakID) = length(rMZ);
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
    mDataCube = zeros( peakNum, hei*wid );
    for i = 1:hei*wid
        if BlkDS.indMap(i) == 1
            curSpec = dataCube(:,i);
            startPos = 1;
            curPeakID = 1;
            for j = 1:length(mIndVec)-1
                if mIndVec(j+1) ~= mIndVec(j)
                    endPos = j;
                    mDataCube(curPeakID,i) = max(dataCube(startPos:endPos, i));
                    startPos = j+1;
                    curPeakID = curPeakID + 1;
                end
            end
            endPos = j+1;
            mDataCube(curPeakID,i) = max(dataCube(startPos:endPos, i));
        end
    end
    mDataCube = reshape( mDataCube, peakNum, hei, wid );
end