function [ dif, dfVec, tVec ] = aggreTestThetaResults( folderName,  seq )
%aggreTestLambdaResults 
cd(folderName);
listing = dir(folderName);
cnt = [];
dfVec = [];
tVec = [];
for i = 1:length(seq)
    fileName = ['testThetaScript_' num2str(i) '_.mat'];

    for j = 3:length(listing)
        if strcmp( listing(j).name, fileName ) == 1
            fprintf( '%s\n', listing(j).name );
            cnt = [cnt i];
            load( fileName );
            dfVec = [dfVec df];
            tVec = [tVec lambda];
        end
    end
end
dif = setdiff(1:length(seq), cnt);

