function [ expNames, outFilePath ] = preprocess2dataCube( inputFileDir, outFileDir)
%preprocess2dataCube
%inputFileDir: The folder containing the preprocessed files
%e.g. each experiment containing three files (if experiment a)
%a_data.csv, a_mz.csv, a_pos.csv
%outFileDir: the folder to output files, file name will be a_dc.mat
fNames = dir( inputFileDir );
expNames = cell( (length(fNames)-2)/4, 1 );
cnt = 1;
for i = 3:(length(fNames)-2)
    res = regexp(fNames(i).name,'(.+)_data.csv','tokens' );
    if ~isempty( res )
        expNames{cnt} = res{1}{1};
        cnt = cnt + 1;
    end
end
outFilePath = cell( length(expNames), 1 );
for z = 1:length( expNames )
    fprintf( ['processing experiment: ', expNames{z}, '\n'] );
    fNamePre = expNames{z};
    dataMatrix = csvread([inputFileDir, '\', fNamePre,'_data.csv']);
    mzAxis = csvread([inputFileDir, '\', fNamePre,'_mz.csv']);
    posInfo = csvread([inputFileDir, '\', fNamePre,'_pos.csv']);
    indMatrix = csvread([inputFileDir, '\', fNamePre,'_indPeak.csv']);
    %whos dataMatrix
    maxXVal = max( posInfo(:, 1) ) + 1;
    maxYVal = max( posInfo(:, 2) ) + 1;
    dataCube = zeros( size(dataMatrix, 1), maxYVal, maxXVal );
    for i = 1:length(posInfo)
        x = posInfo(i, 1);
        y = posInfo(i, 2);
        dataCube(:, y, x) = dataMatrix(:, i);
    end
    dataCube = round(dataCube);
    outFilePath{z} = [outFileDir, '\', fNamePre, '_dc.mat'];
    IMSData.dataCube = dataCube; IMSData.mzAxis = mzAxis; IMSData.posInfo = posInfo; IMSData.indMatrix = indMatrix;
    save( outFilePath{z}, 'IMSData', '-v7.3' );
end
end

