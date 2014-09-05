function [ gridVec ] = genHypParamGrid( upperVal, gridNum )
%genHypParamGrid generate hyper parameter grid for different hyper
%parameters learning
%input
%upperVal a vector of dim 2 or 3, m is # hyper parameters. This vector contains the
%upper bound of hyper parameters
%gridNum a scalr, the number of ticks
%output
%depends on dim, if dim = 3, then gridVec has three fields: fVec, sVec,
%tVec.
dNum = length( upperVal );
if dNum == 3
    stepVec = zeros( 3, 1 );
    stepVec(1) = upperVal(1)/gridNum; stepVec(2) = upperVal(2)/gridNum; stepVec(3) = upperVal(3)/gridNum;
    [fVec, sVec, tVec] = meshgrid( stepVec(1):stepVec(1):upperVal(1)...
        , stepVec(2):stepVec(2):upperVal(2), stepVec(3):stepVec(3):upperVal(3) );
    gridVec.fVec = fVec; gridVec.sVec = sVec; gridVec.tVec = tVec;
elseif dNum == 2
    stepVec = zeros( 2, 1 );
    stepVec(1) = upperVal(1)/gridNum; stepVec(2) = upperVal(2)/gridNum;
    [fVec, sVec] = meshgrid( stepVec(1):stepVec(1):upperVal(1)...
        , stepVec(2):stepVec(2):upperVal(2) );
    gridVec.fVec = fVec; gridVec.sVec = sVec;
end
end

