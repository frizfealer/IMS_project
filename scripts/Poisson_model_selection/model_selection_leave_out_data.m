function [ WRes1, WRes2 ] = model_selection_leave_out_data( dataCube, expRecReal, DTemplate, SpeciesM, DIonName, mzAxis, type, outName, refRegionName  )
%model_selection_leave_out_data model, select model of origianl
%dictionary element d and elements that generated from d
%dataCube: experiment data
%expRecReal: experiment result
%DTemplate: dictionary template
%SpeciesM: species M for each element of DTemplate
%DIonName: Ion name for each element of DTemplate
%mzAxis: mzAxis for the experiment
%type: 'facilitate' or 'hinder'
%outName: resulting file path
%refRegionName: reference region of interaction colonies

%% default constants:
%the intensity of dictionary entries to be ignore
IGNORE_INT = 1e-2;
%the number of element need to be tested
TEST_DIC_ELE_NUM = 30;

BlkDS = conBLKDS( dataCube );
[~, ~, sigIdx] = analyzeResults( expRecReal, DTemplate, BlkDS, SpeciesM, DIonName, mzAxis, [], 892, '~/', 0, type, refRegionName );
D = expRecReal.outD;
[ aMatrix ] = genAMatrix( BlkDS, 0.1, 'byRow' );
[~, hei, wid] = size( dataCube );
nLen = hei*wid;

%% set z0
z0 = log(dataCube);
z0(z0==-inf)=0;
%if not in traing set, we shoud not initialize z0 according to it.
%take the average values
tmp = sum( z0(:, :), 2 ) / ( length( find( BlkDS.indMap == 1 ) ) );
for i = 1:hei*wid
    %if is in data area, but not in training set, 
    %set it to average spectrum
    if BlkDS.indMap(i) && ~aMatrix(i)
        z0(:, i) = tmp;
    end
end
%% set scale factor
scaleFactor =  1 / ( max( dataCube(:) ) / max( z0(:) ) ) * 100;

%% set phi theta
phi = expRecReal.phi;
theta = expRecReal.theta;

%% test for two hypothesis for each dictionary elements
%The first is one dictionary elements
%another is n dictionary elements (assume a dictionary element with n
%entries > IGNORE_INT), and each elements only has one entry.
for i = 1:TEST_DIC_ELE_NUM
    cMZ = find( D(:, sigIdx(i)) > IGNORE_INT );
    if length(cMZ) == 1
        continue;
    end
    
    %% first model
    firstD = D(cMZ, sigIdx(i));
    W0 = zeros( 1, nLen );
    W = zeros( 1, nLen );
    z1 = zeros( size( W ) );
    [res] = updateW_ADMM...
        ( dataCube(cMZ,:, :), firstD, W, W0, z0(cMZ, :, :), z1, aMatrix, BlkDS, 100, 0, phi, theta, scaleFactor, [] );
    WRes1(i) = res;
    
    %% second model
    secondD = zeros( length(cMZ) );
    for j = 1:length(cMZ)
        secondD(j, j) = D(cMZ(j), sigIdx(i));
    end
    W = zeros( length(cMZ), nLen );
    W0 = zeros( 1, nLen );
    z1 = zeros( size( W ) );
    [res] = updateW_ADMM...
        ( dataCube(cMZ,:, :), secondD, W, W0, z0(cMZ, :, :), z1, aMatrix, BlkDS, 100, 0, phi, theta, scaleFactor, [] );
    WRes2(i) = res;

    save( outName, 'WRes1', 'WRes2' );
end



end

