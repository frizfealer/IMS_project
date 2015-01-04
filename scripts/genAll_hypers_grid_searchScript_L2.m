function [ ] = genAll_hypers_grid_searchScript_L2( folderName )
%genTestLambdaScript generate TestLambdaScript for each seq
mkdir(folderName);
cd(folderName);
seq = 1:125;
for i = 1:length(seq)
    fileName = ['all_hypers_grid_search_L2Script_' num2str(i) '_.m'];
    fileID = fopen(fileName, 'w+');
    fprintf( fileID, 'addpath( genpath( ''~/IMS_project-master'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/glmnet_matlab'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/Ipopt_3118'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/SLEP_package_4.1'' ) );\n' );
    fprintf( fileID, 'addpath( genpath( ''~/minFunc_2012'' ) );\n' );
    fprintf( fileID, ['exp_100831_348_136_12_40hr_0_1XLB_LN_all_hypers_grid_search_L2( ' num2str(i) ' );\n'] );
    fclose(fileID);
end


end

