function [ dataCube, mzAxis ] = preprocess2dataCube( dataFilePath, mzAFilePath, pDataFilePath, outFilePath )
%preprocess2dataCube
dataMatrix = csvread(dataFilePath);
mzAxis = csvread(mzAFilePath);
posInfo = csvread(pDataFilePath);
whos dataMatrix
maxXVal = max( posInfo(:, 1) ) + 1;
maxYVal = max( posInfo(:, 2) ) + 1;
dataCube = zeros( size(dataMatrix, 1), maxYVal, maxXVal );
for i = 1:length(posInfo)
    x = posInfo(i, 1);
    y = posInfo(i, 2);
    dataCube(:, y, x) = dataMatrix(:, i);
end
if ~isempty( outFilePath )
    save( outFilePath, 'dataCube', 'mzAxis' );
end
end

