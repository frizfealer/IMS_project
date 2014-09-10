function [ aMatrix ] = genAMatrix( BlkDS, percent, methods)
%genAMatrix generate alpha Matrix for leave out testing
%BlkDS, block data structure from conBLKDS
%percent, the percentage of data in testing
%byRow, is == 1, choose testing samples by row, from the upper-left corner
%output:
%aMatrx the same size as gW, with 1 means trainning and 0 means two
%situations: one is that location has no sample, another is that location
%is for testing. To distinguish two cases, use aMatrix==0&BlkDS.indMap==1
%to show the testing cases.
[hei, wid] = size( BlkDS.indMap );
aMatrix = zeros( size( BlkDS.indMap ) );
if strcmp( methods, 'byRow' ) == 1
    mFlag = 1;
elseif strcmp( methods, 'byCol' ) == 1
    mFlag = 0;
elseif strcmp( methods, 'grid' ) == 1
    mFlag = 2;
end
GRIDLEN = 5;
for i = 1:BlkDS.blkNum
    loc = BlkDS.B2GMap{i};
    aMatrix(loc) = 1;
    [I,J] = ind2sub( [hei wid], loc );
    [I, idx] = sort(I);
    J = J(idx);
    aNum = round( length(loc)*percent );
    if mFlag == 1
        aMatrix(I(1:aNum), J(1:aNum)) = 0; 
    elseif mFlag == 0
        aMatrix(loc(1:aNum)) = 0; 
    elseif mFlag == 2
        
    end
end
aMatrix = aMatrix == 1;

end

