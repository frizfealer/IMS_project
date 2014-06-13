function [ tBlkDS ] = conBLKDS( dataCube )
%% conDSBlkIdx construct Data Structure of block information
%Input:
%bImg, a binary image
%Output: a data structure conatins the following fields:
%blkNum: # blocks
%B2WMap: Block id to Grid coordinates mapping
%G2BMap: Grid coordinate to Block ids mapping
%indMap: indicator for all location in the grid, with 1 means computation
%needed.

bImg = sum( dataCube, 1 );
[~, HEI, WID] = size( dataCube );
bImg = reshape( bImg, HEI, WID ); bImg( bImg > 0 ) = 1;
%bwboundaries, return [y,x]
[B,L] = bwboundaries( bImg );

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

