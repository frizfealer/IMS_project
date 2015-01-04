function [ ] = genAll_hypers_grid_searchScript_LNR( finFilesFolder, folderName )
%genTestLambdaScript generate TestLambdaScript for each seq
mkdir(folderName);
cd(folderName);
listing = dir(finFilesFolder);
finNum = [];
for i = 3:length(listing)
    cNum = regexp( listing(i).name, ...
        'all_hypers_grid_search_resultsexp_100831_348_136_12_40hr_0_1XLB_LN_res_(\d+)_.mat', 'tokens');
    tmp = cNum{1}; 
    tmp = str2double(tmp);
    finNum = [finNum  tmp];
end
unfinNum = setdiff( 1:125, finNum );
for i = 1:length(unfinNum)
    fileName = ['all_hypers_grid_searchScript_LNR_' num2str(unfinNum(i)) '_.m'];
    fileID = fopen(fileName, 'w+');
    fprintf( fileID, 'addpath( genpath( ''~/IMS_project-master'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/glmnet_matlab'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/Ipopt_3118'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/SLEP_package_4.1'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/minFunc_2012'' ) );\n' );
    fprintf( fileID, ['exp_100831_348_136_12_40hr_0_1XLB_LN_all_hypers_grid_search_LNR( ' num2str(unfinNum(i)) ' );\n'] );
    fclose(fileID);
end

end

