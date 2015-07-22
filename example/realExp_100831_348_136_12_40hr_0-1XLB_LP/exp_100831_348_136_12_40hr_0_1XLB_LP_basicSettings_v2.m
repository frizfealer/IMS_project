function [ expRec, gridVec ] = exp_100831_348_136_12_40hr_0_1XLB_LP_basicSettings_v2( INPUT_FOLDER_PATH, OUTPUT_FILE_PATH)
%experiment scripts update @ 2015/07/21
%PROJECT_FOLDER_PATH: you should change to your own folder of the codes.
PROJECT_FOLDER_PATH = 'D:\Users\YeuChern\GitHub\IMS_project';
%INPUT_FOLDER_PATH = 'D:\IMS_DATA\100831_348_136_12,40hr_0-1XLB_LP\outCSV';
%OUTPUT_FILE_PATH = '\example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12_40hr_0_1XLB_LP_vars.mat';
%% convert four csv files into a datacube and save it to the experiment folder
[ ~, ~ ] = preprocess2dataCube( INPUT_FOLDER_PATH, ...
    [PROJECT_FOLDER_PATH,'\example\realExp_100831_348_136_12_40hr_0-1XLB_LP']);

%% load data, including dataCube and mzAxis
load( [PROJECT_FOLDER_PATH,'\example\realExp_100831_348_136_12_40hr_0-1XLB_LP\100831_348_136_12,40hr_0-1XLB_LP.mzML_dc.mat'] );

%% generate BlkDS, a data structure for bacteria community location information
%BlkDS has the 4 fields
%blkNum: the bacteria community number 
%B2GMap: a mapping from bacteria community to grid
%G2BMap: a mapping from grid to bacteria community
%indMap: a logical matrix, as the same size of the grid, with 1 means the
%grids having signals, 0 means empty
IMSD = IMSData; clear IMSData;
IMSD.BlkDS = conBLKDS( IMSD.dataCube );

%% bin dataCube
[ IMSD.dataCube, mMZAxis, IMSD.mappingFunc, ~, ~, ~ ] = binningDataCube( IMSD.dataCube, IMSD.mzAxis, IMSD.BlkDS, []);
IMSD.oMZAxis = IMSD.mzAxis; IMSD.mzAxis = mMZAxis; clear mMZAxis;

%% retrieve data dimension
[SLEN, IHEIGHT, IWIDTH] = size( IMSD.dataCube );

%% generate DTemplate
IonTableFilePathPos = 'D:\Users\YeuChern\GitHub\IMS_project\example\molecule_profile_pos_v2.csv';
%[ DTemplate2, DIonName2, speciesM2 ] = genDTemplate_v3( IMSD.mzAxis, IonTableFilePathPos, 5e-4 );
%[ DTemplate2, DIonName2, SpeciesM2 ] = improveDTemplate( 'positive', DTemplate2, DIonName2, speciesM2 );
[ DTemplate, DIonName, speciesM ] = genDTemplate_v4( IonTableFilePathPos, IMSD, 5e-4 );
[ DTemplate, DIonName ] = addingIsotopeDTemplate( IMSD.mzAxis, 5e-4, DTemplate, DIonName, 2, [] );

%% generate leave-out pattern among grid cells 
[ aMatrix ] = leaveoutCBPatternsData_v2( IMSD.BlkDS, 10 );

%% saves
save( OUTPUT_FILE_PATH, 'aMatrix', 'DTemplate', 'DIonName', 'speciesM', 'SLEN', 'IHEIGHT', 'IWIDTH', 'IMSD' );

end


