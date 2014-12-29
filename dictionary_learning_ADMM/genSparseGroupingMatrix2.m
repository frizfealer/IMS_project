function [ R ] = genSparseGroupingMatrix2( inMat, m, offsetFlag )
%genSparseGroupingMatrix2 generate sparse grouping matrix for 2D fused lasso
%term. The input comes from sparseGroupingMatrix from each colonies' areas.
%The output is for W vectorized (W(:))
%% for the first type of relation: horizontal fused terms
relNum = size(inMat, 1);
yAxis = zeros( relNum*2*m, 1 );
xAxis = zeros( relNum*2*m, 1 );
valVec = zeros( relNum*2*m, 1);
if offsetFlag == 0
    for i = 1:relNum
        nonZIdx = find( abs( inMat(i, :) ) == 1 );
        yVal = repmat( (1+(i-1)*m):(i*m), 2, 1); 
        yVal = yVal(:);
        yAxis( (1+(i-1)*2*m):(i*2*m) ) = yVal;
        xVal = [(1+(nonZIdx(1)-1)*m):(nonZIdx(1)*m); (1+(nonZIdx(2)-1)*m):(nonZIdx(2)*m)];
        xVal = xVal(:);
        xAxis( (1+(i-1)*2*m):(i*2*m) ) = xVal;
        valVal = repmat( inMat(i, nonZIdx ), 1, m );
        valVal = valVal(:);
        valVec( (1+(i-1)*2*m):(i*2*m) ) = valVal;
    %     yAxis( (1+(i-1)*2*m):(i*2*m) ) = i;
    %     xAxis( (1+(i-1)*2*m):(i*2*m-m) ) = (1+(nonZIdx(1)-1)*m):(nonZIdx(1)*m);
    %     xAxis( (i*2*m-m+1):(i*2*m) ) = (1+(nonZIdx(2)-1)*m):(nonZIdx(2)*m);
    %     valVec( (1+(i-1)*2*m):(i*2*m-m) ) = repmat( inMat(i, nonZIdx(1) ), m, 1 );
    %     valVec( (i*2*m-m+1):(i*2*m) ) = repmat( inMat(i, nonZIdx(2) ), m, 1 );
    end
    R = sparse( yAxis, xAxis, valVec, relNum*m, size(inMat, 2)*m );
else
    m = m + 1;
    for i = 1:relNum
        nonZIdx = find( abs( inMat(i, :) ) == 1 );
        yVal = repmat( (1+(i-1)*m):(i*m), 2, 1);
        yVal = yVal(:);
        yAxis( (1+(i-1)*2*m):(i*2*m) ) = yVal;
        xVal = [(1+(nonZIdx(1)-1)*m):(nonZIdx(1)*m); (1+(nonZIdx(2)-1)*m):(nonZIdx(2)*m)];
        xVal = xVal(:);
        xAxis( (1+(i-1)*2*m):(i*2*m) ) = xVal;
        valVal = repmat( inMat(i, nonZIdx ), 1, m-1 );
        valVal = [valVal zeros(1, 2)];
        valVal = valVal(:);
        valVec( (1+(i-1)*2*m):(i*2*m) ) = valVal;
        %     yAxis( (1+(i-1)*2*m):(i*2*m) ) = i;
        %     xAxis( (1+(i-1)*2*m):(i*2*m-m) ) = (1+(nonZIdx(1)-1)*m):(nonZIdx(1)*m);
        %     xAxis( (i*2*m-m+1):(i*2*m) ) = (1+(nonZIdx(2)-1)*m):(nonZIdx(2)*m);
        %     valVec( (1+(i-1)*2*m):(i*2*m-m) ) = repmat( inMat(i, nonZIdx(1) ), m, 1 );
        %     valVec( (i*2*m-m+1):(i*2*m) ) = repmat( inMat(i, nonZIdx(2) ), m, 1 );
    end
    R = sparse( yAxis, xAxis, valVec, relNum*m, size(inMat, 2)*m );
end
end
