function [ ] = gen_hypers_grid_search_script( folderName, fileNamePrefix, exeFuncName, seqNum )
%--------------------------------------------------------------------------
% gen_hypers_grid_search_script: a script to generate all scripts to do
% the hyper parameters grid search in a cluter.
%--------------------------------------------------------------------------
% DESCRIPTION:
%   This is a script to generate a lot of scripts with different hyper
%   parameters. The reason to do this is this can run all hyper-parameters
%   searching at a time. This is done in our cluster (killdevil.unc.edu)
% INPUT ARGUMENTS:
%   folderName, the output folder path name.
%   fileNamePrefix, the output file name prefix, and there will be a
%   "_number" adding to the prefix.
%   exeFuncName, the execution function name.
%   seqNum, the hyper parameters index.


oPath = pwd;
mkdir(folderName);
cd(folderName);
seq = 1:seqNum;
for i = 1:length(seq)
    fileName = [fileNamePrefix '_' num2str(i) '.m'];
    fileID = fopen(fileName, 'w+');
    fprintf( fileID, 'addpath( genpath( ''~/IMS_project-master'' ) );\n' );

%     fprintf( fileID, ['exp_100831_348_136_12_40hr_0_1XLB_LN_two_hypers_grid_search( ' num2str(i) ' );\n'] );
    fprintf( fileID, [exeFuncName '( ' num2str(i) ' );\n'] );
    fclose(fileID);
end
cd( oPath );

end