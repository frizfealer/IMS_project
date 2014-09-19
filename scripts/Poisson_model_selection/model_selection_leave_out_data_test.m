function [ resMat ] = model_selection_leave_out_data_test( dataCube, expRecReal, DTemplate, SpeciesM, DIonName, mzAxis, type, refRegionName, inModelsPath  )
%model_selection_leave_out_data model_test, tset model of origianl
%dictionary element d and elements that generated from d on leave-out-data
%dataCube: experiment data
%expRecReal: experiment result
%DTemplate: dictionary template
%SpeciesM: species M for each element of DTemplate
%DIonName: Ion name for each element of DTemplate
%mzAxis: mzAxis for the experiment
%type: 'facilitate' or 'hinder'
%outName: resulting file path
%refRegionName: reference region of interaction colonies
%inModelsPath: the *.mat file generated from model_selection_leave_out_data

load( inModelsPath );
BlkDS = conBLKDS( dataCube );
[~, ~, sigIdx] = analyzeResults( expRecReal, DTemplate, BlkDS, SpeciesM, DIonName, mzAxis, [], 892, '~/', 0, type, refRegionName );
D = expRecReal.outD;
IGNORE_INT = 1e-2;
[ aMatrix ] = genAMatrix( BlkDS, 0.1, 'byRow');
[~, hei, wid] = size(dataCube);

%% compute scaleFactor
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

%% compute test results
DIFF_VAL = 5*1e-2;

testArea = aMatrix==0&BlkDS.indMap==1;
testLL = zeros( 30, 4 );
testMean = zeros( 30, 2 );
testStd = zeros( 30, 2 );
locMat = zeros(30, 10 );
resMat = zeros( 1, 7 );
for i = 1:length(WRes1)
    fprintf('%d\t', i);
    cMZ = find( D(:, sigIdx(i)) > IGNORE_INT );
    if length(cMZ) == 1
        continue;
    end
    
    %% first model
    firstD = D(cMZ, sigIdx(i));
    W = WRes1(i).W;
    
    [centers] = determineBCen( W(1, BlkDS.indMap), DIFF_VAL );
    [dW] = discretizeW( W, centers, BlkDS.indMap );
    W0 = WRes1(i).W0;
    %[ val1, meanVal1, stdVal1, ] = LP_DL_Poiss( testArea, dataCube(cMZ, :), W, W0, firstD, 0, 0, 0, scaleFactor, [], 1 );
    LLVal1 = easy_PoissonLL( testArea, dataCube(cMZ, :), firstD, W, W0 );
    [ df1 ] = findNumConCom( dW );
    %% second model
    secondD = zeros( length(cMZ) );
    for j = 1:length(cMZ)
        secondD(j, j) = D(cMZ(j), sigIdx(i));
    end
    W = WRes2(i).W;
    dW = W;
    df2 = zeros( length(cMZ), 1 );
    for j = 1:length(cMZ)
        [centers] = determineBCen( W(j, BlkDS.indMap), DIFF_VAL );
        [dW(j, :)] = discretizeW( W(j, :), centers, BlkDS.indMap );
        df2(j) = findNumConCom(dW(j,:));
    end
    W0 = WRes2(i).W0;
    %[ val2, meanVal2, stdVal2 ] = LP_DL_Poiss( testArea, dataCube(cMZ, :), W, W0, secondD, 0, 0, 0, scaleFactor, [], 1 );
    LLVal2 = easy_PoissonLL( testArea, dataCube(cMZ, :), secondD, W, W0 );
    [h, pVal] = lratiotest( LLVal2,LLVal1,sum(df2) - df1,0.01);    
    fprintf( 'LLVal2 = %g, LLVal1 = %g, df2 = %d, df1 = %d, h = %d, pVal = %g\n',LLVal2, LLVal1, sum(df2), df1,  h, pVal );
    resMat(i,:) = [i, LLVal2, LLVal1, sum(df2), df1, h, pVal];
    %testLL(i, :) = [val1 val2, LLVal1, LLVal2];
    %testMean(i, :) = [meanVal1 meanVal2];
    %testStd(i, :) = [stdVal1 stdVal2];
end

end

function [centers] = determineBCen( W, diff_val )
binNum = 100;
while 1
    [~,centers] = hist( W, binNum );
    temp = diff(centers);
    if temp(1) <= diff_val
        break;
    end
    binNum = binNum + 5;
end
end

function [dW] = discretizeW( W, centers, indMap )
        %discretize
    dW = zeros( size(W) );
    for j = 1:length(W(:))
        if indMap(j) == 1
            [~,idx] = min( abs( centers-W(j) ) );
            dW(j) = centers(idx);
        end
    end
end

function [ val ] = easy_PoissonLL( aMatrix, Y, D, W, W0)
Y = Y(:, :);
D = D(:, :);
W = W(:, :);
[sLen, ~] = size( Y );

preY = D* W + repmat( W0(:)', sLen, 1 ); %preY [s ,w*h]
lFY = logfactorial_e( Y(:, aMatrix), 1e8 );
resMat = Y(:, aMatrix).*preY(:, aMatrix) - exp( preY(:, aMatrix) ) - lFY;
val = sum(sum(resMat));
end

function [ val ] = findNumConCom( dW )
    WVal = unique(dW);
    val = 0;
    for i = 1:length(WVal)
        if WVal(i) == 0
            continue;
        end
        BW = full(dW);
        BW( BW ~= WVal(i) ) = 0;
        B = bwboundaries( BW, 8, 'noholes' );
        val = val + length(B);
    end
end