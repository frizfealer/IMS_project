function [ output_args ] = model_selection_leave_out_data( dataCube, expRecReal, pDTemplate, pSpeciesM, pDIonName, mzAxis  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
BlkDS = conBLKDS( dataCube );
[ sigDict, sigM, sigIdx] = analyzeResults( expRecReal, pDTemplate, BlkDS, pSpeciesM, pDIonName, mzAxis, [], 892, '~/', 0 );
D = expRecReal.outD;
IGNORE_INT = 1e-2;
[ aMatrix ] = genAMatrix( BlkDS, 0.1, 'byRow');
[sLen, hei, wid] = size(dataCube);
nLen = hei*wid;

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
scaleFactor =  1 / ( max( dataCube(:) ) / max( z0(:) ) ) * 100;

phi = expRecReal.phi;
theta = expRecReal.theta;

%% test for two hypothesis for each dictionary elements
%one is one dictionary elements
%another is n dictionary elements (assume a dictionary element with n
%entries > IGNORE_INT), nd each elements only has one entry.
TEST_DIC_ELE_NUM = 30;
for i = 14:TEST_DIC_ELE_NUM
    cMZ = find( D(:, sigIdx(i)) > IGNORE_INT );
    if length(cMZ) == 1
        continue;
    end
    firstD = D(cMZ, sigIdx(i));
    W0 = zeros( 1, nLen );
    W = zeros( 1, nLen );
    z1 = zeros( size( W ) );
    [res] = updateW_ADMM...
        ( dataCube(cMZ,:, :), firstD, W, W0, z0(cMZ, :, :), z1, aMatrix, BlkDS, 100, 0, phi, theta, scaleFactor, [] );
    WRes1(i) = res;
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

    save( 'res_differentModels_continue.mat', 'WRes1', 'WRes2' );
end



end

