function [ tBlkDS ] = conBLKDS( dataCube )
%--------------------------------------------------------------------------
% conDSBlk construct Data Structure of block information
%--------------------------------------------------------------------------
% DESCRIPTION:
%   Since there is only a subset of grids in the data cube having data, we
%   only want to train our model on this subset of grids.
%
% INPUT ARGUMENTS:
%   dataCube, size of [s w h], that is #(m/z)*width*height
% OUTPUT ARGUMENTS:
%   tBlkDS a data structrue having fields:
%       blkNum: # communities
%       B2WMap: a mapping from a block ID to grid coordinates
%       G2BMap: a mapping from a grid coordinate to a block ID, with empty
%       means no data in the coordinate.
%       indMap: indicator for all grids with one means having signals
%       (computation needed) and zero means no signals.

bImg = sum( dataCube, 1 );
[~, HEI, WID] = size( dataCube );
bImg = reshape( bImg, HEI, WID ); bImg( bImg > 0 ) = 1;
%bwboundaries, return [y,x]
[B,L] = bwboundaries( bImg, 8, 'noholes' );

blkNum = length(B);

B2GMap = cell( blkNum, 1 );
for j = 1:blkNum
    %find(L == j)
    %linearInd = sub2ind( [HEI WID], L{1}(:,1), L{1}(:,2));
    B2GMap{j} = find(L == j);
end
%% construct W2BMap
G2BMap = cell( HEI*WID, 1 );
for i = 1:HEI*WID
    G2BMap{i} = [];
end
for i = 1:blkNum
    cBIdx = B2GMap{i};
    for j = cBIdx'
        G2BMap{j} = [G2BMap{j} i]; 
    end
end

tBlkDS.blkNum = blkNum;
tBlkDS.B2GMap = B2GMap;
tBlkDS.G2BMap = G2BMap;
tBlkDS.indMap = bImg == 1;

end

