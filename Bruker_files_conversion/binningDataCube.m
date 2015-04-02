function [ mDataCube, mMZAxis, mappingFunc, mMZAxis2, mappingFunc2, tolVec ] = binningDataCube( dataCube, mzAxis, BlkDS, posBinMZ)
%--------------------------------------------------------------------------
% binningDataCube: binning data Cube with choosing the conesus peak
% positions among samples and wiht TOL = 5e-4 (tolerance = 500 ppm)
%--------------------------------------------------------------------------
% DESCRIPTION:
%   Since there are peak shifting across samples grids, we do binning on
%   the data cube.
%
% INPUT ARGUMENTS:
%   dataCube, size of [s w h], that is #(m/z)*width*height
%   mzAxis, size of [s 1], a vector of m/z values of the input data.
%   BlkDS, block structure from conBLKDS
%   posBinMZ, deprecated now.
% OUTPUT ARGUMENTS:
%   mDataCube, the data cube after binning
%   mMZAxis, the mzAxis after binning
%   mappingFunc, a mapping from the a set of indices in mzAxis to their binned m/z value in
%   mMZAxis

TOL = 5e-4;

[~, hei, wid] = size(dataCube);
ins = dataCube(:, BlkDS.indMap == 1 );
ins = sum( ins, 2 );
[~,dOrd] = sort(ins, 'descend');
dMZ = mzAxis(dOrd);
mMZAxis = zeros( length(dMZ), 1 );
mIndVec = zeros( length(dMZ), 1 );
curPeakID = 1;
tolVec = [];
while ~isempty(dMZ)
    %         fprintf( '%d\n', curPeakID );
    %         idx = ismember( mzAxis, dMZ ) ;
    %         tmp = ins(idx);
    %         if max(tmp) ~= ins( mzAxis == dMZ(1) )
    %             keyboard();
    %         end
    cTOL = TOL;
    lMZ = dMZ(1);
    if ~isempty( posBinMZ )
        if ismember(lMZ, posBinMZ)
            cTOL = 1e-4;
        end
    end
    tolVec = [tolVec; cTOL];
    rMZ = dMZ( dMZ >= lMZ - lMZ*cTOL );
    rMZ = rMZ( rMZ <= lMZ + lMZ*cTOL );
    mIndVec( ismember( mzAxis, rMZ )==1 ) = curPeakID;
    mMZAxis(curPeakID) = lMZ;
    %iii(curPeakID) = length(rMZ);
    curPeakID = curPeakID + 1;
    rDOrd = ismember( dMZ, rMZ );
    %         fprintf( '%d, %d, %d\n', length(dMZ), length(rMZ), length(rDOrd) );
    %         if length(rMZ) ~= length(rDOrd)
    %             keyboard();
    %         end
    dMZ(rDOrd) = [];
end

peakNum = max(mIndVec);
tmp = mMZAxis;
mMZAxis = tmp(1:peakNum);
[mMZAxis, idx] = sort( mMZAxis );
tolVec = tolVec(idx);
mMZAxis2 = zeros( size( mMZAxis ) );
mDataCube = zeros( peakNum, hei*wid );
%     binMap = zeros( size( mDataCube ) );
mappingFunc = zeros( size(mIndVec) );
mappingFunc2 = zeros( size(mIndVec) );
cnt = 0;
for i = 1:hei*wid
    if BlkDS.indMap(i) == 1
        cnt = cnt + 1;
        startPos = 1;
        curPeakID = 1;
        for j = 1:length(mIndVec)-1
            if mIndVec(j+1) ~= mIndVec(j)
                endPos = j;
                mDataCube(curPeakID,i) = max(dataCube(startPos:endPos, i));
                mappingFunc(startPos:endPos) = mMZAxis(curPeakID);
                %[~, idx] = max(dataCube(startPos:endPos, i));
                %                     binMap(curPeakID, i) = max( mzAxis(idx+startPos-1)- mzAxis(startPos), mzAxis(endPos)-mzAxis(idx+startPos-1) );
                mMZAxis2(curPeakID) = mean( mzAxis(startPos:endPos) );
                mappingFunc2(startPos:endPos) = mMZAxis2(curPeakID);
                %                     binMap(curPeakID, i) = max( mMZAxis2(curPeakID)- mzAxis(startPos), mzAxis(endPos)-mMZAxis2(curPeakID) );
                startPos = j+1;
                curPeakID = curPeakID + 1;
            end
        end
        endPos = j+1;
        mDataCube(curPeakID,i) = max(dataCube(startPos:endPos, i));
        mMZAxis2(curPeakID) =mean( mzAxis(startPos:endPos) );
        mappingFunc2(startPos:endPos) = mMZAxis2(curPeakID);
        mappingFunc(startPos:endPos) = mMZAxis(curPeakID);
    end
end
mDataCube = reshape( mDataCube, peakNum, hei, wid );
end