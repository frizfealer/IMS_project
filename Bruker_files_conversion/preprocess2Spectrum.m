function [] = preprocess2Spectrum( inputFilePath, outFileDir )
%preprocess2Spectrum for one spectrum file
res = regexp(inputFilePath,'(.+)\\(.+)_data.csv','tokens' );
fNamePre = res{1}{2};
inputFileDir = res{1}{1};
dataVec = csvread([inputFileDir, '\', fNamePre,'_data.csv']);
mzAxis = csvread([inputFileDir, '\', fNamePre,'_mz.csv']);
outFilePath = [outFileDir, '\', fNamePre, '_dc.mat'];
save( outFilePath, 'dataVec', 'mzAxis', '-v7.3' );
end

