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

testArea = aMatrix==0&BlkDS.indMap==1;
testLL = zeros( 30, 2 );
testMean = zeros( 30, 2 );
testStd = zeros( 30, 2 );
for i = 1:30
    fprintf('%d\n', i);
    cMZ = find( D(:, sigIdx(i)) > IGNORE_INT );
    if length(cMZ) == 1
        continue;
    end
    firstD = D(cMZ, sigIdx(i));
    W = WRes1_c(i).W;
    W0 = WRes1_c(i).W0;
    [ val1, meanVal1, stdVal1, ] = LP_DL_Poiss( testArea, dataCube(cMZ, :), W, W0, firstD, 0, 0, 0, scaleFactor, [], 1 );
    secondD = zeros( length(cMZ) );
    for j = 1:length(cMZ)
        secondD(j, j) = D(cMZ(j), sigIdx(i));
    end
    W = WRes2_c(i).W;
    W0 = WRes2_c(i).W0;
    [ val2, meanVal2, stdVal2 ] = LP_DL_Poiss( testArea, dataCube(cMZ, :), W, W0, secondD, 0, 0, 0, scaleFactor, [], 1 );
    testLL(i, :) = [val1 val2];
    testMean(i, :) = [meanVal1 meanVal2];
    testStd(i, :) = [stdVal1 stdVal2];
end

for i = 1:30
outString1{i}=['Model 1: ', num2str(testMean(i,1)), ' +/- ', num2str(testStd(i,1))];
outString2{i}=['Model 2: ', num2str(testMean(i,2)), ' +/- ', num2str(testStd(i,2))];
end

fd = fopen( 'diffModel.txt', 'w+');
for i = 1:30
fprintf( fd, '%s\t%s\n', outString1{i}, outString2{i} );
end
fclose( fd );