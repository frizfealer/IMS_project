function [ dataCube, mzAxis ] = preprocess2dataCube( inputFilePath, outFilePath )
%preprocess2dataCube
%inputFilePath: The folder containing three input files e.g. D:
%outFilePath: the folder and the file name (.mat) of output file e.g.
%D:\out.mat
dataMatrix = csvread([inputFilePath,'\dataMatrix.csv']);
mzAxis = csvread([inputFilePath,'\mzVec.csv']);
posInfo = csvread([inputFilePath,'\posMat.csv']);
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

