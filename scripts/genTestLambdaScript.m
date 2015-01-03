function [ ] = genTestLambdaScript( folderName, seq )
%genTestLambdaScript generate TestLambdaScript for each seq
mkdir(folderName);
cd(folderName);
for i = 1:length(seq)
    fileName = ['testLambdaScript_' num2str(i) '_.m'];
    fileID = fopen(fileName, 'w+');
    fprintf( fileID, 'addpath( ''/csbiohome01/ycharn '');\n');
    fprintf( fileID, 'startup\n' );
    fprintf( fileID, 'load( ''/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_input_20141231.mat'' );\n');
    fprintf( fileID, 'iD = initD( dataCube, nDTemplate, ''NNMF'', nDIonName);\n' );
    fprintf( fileID, ['[df, lambda] = testLambdaMax( dataCube, iD, 1, 1, 1e-1, 1,' num2str(seq(i)) ', 5e-2 );\n'] );
    outFileName = [fileName(1:(end-1)) 'mat'];
    fprintf( fileID, 'mkdir(''results'');\n' );
    fprintf( fileID, ['save( ''results/' outFileName ''', ''df'', ''lambda'' )'] );
    fclose(fileID);
end


end

