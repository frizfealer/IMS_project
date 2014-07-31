function [ rD ] = updateD_v5( inY, outW, outW0, D_init, DTemplate, aMatrix, HesOpt, phi, scaleFac, itNum )
%updateD_v5 update the whole dictionary using fmincon without ADMM
%aMatrix, a indicator matrix, with 1 means using in trainning and 0 means
%using in testing
%HesOpt, a flag of using Hessian or not, 1 means using.
%phi, the sparsity weighting parameter
%scaleFac, scalar, to make phi more controlable
[sLen, mLen] = size(D_init);
%% option of fmincon section
options = optimoptions('fmincon');
options.MaxFunEvals= 1e6;
options.MaxIter = itNum;
options.TolFun=1e-16;
options.TolX = 1e-3;
options.TolCon = 1e-16;
% options.TolCon=1e-8;
options.Algorithm='interior-point';
options.Display='iter';
% options.AlwaysHonorConstraints = 'none';
options.UseParallel='always';
options.GradObj='on';
if HesOpt == 1
    options.Hessian = 'user-supplied';
else
    options.Hessian = {'lbfgs',1500};
% options.Hessian = 'fin-diff-grads'; options.SubproblemAlgorithm = 'cg';
end
options.GradConstr = 'on';
% options.PlotFcns = @optimplotfval;
options.TolProjCG = 1e-8;
options.TolProjCGAbs = 1e-16;
options.ScaleProblem='obj-and-constr';
options.InitBarrierParam = max(inY(:));
% options.DerivativeCheck = 'on';

%% initialize variables
aIdx = aMatrix(:)==1;
Y = inY(:, aIdx);
W = outW(:, aIdx);
staticTerm = repmat( outW0(aIdx)', sLen, 1 );
rD = sparse( D_init );

%% generating metadata for nonZeroPosition
nonZPos = find(DTemplate == 1 );
if HesOpt == 1
    fprintf( 'Generating metaData of non-zero positions in D...\n' );
end
%generate nonZPos (one-dimension indicator), map nonzero element -> element
%in the dictionary. nonZPosY (two-dimension indicator, y-axis), map nonzero
%element -> element's Y position in the dictionary
%nonZPosX (two-dimension indicator, x-axis) map nonzero
%element -> element's X position in the dictionary
%nonZPos = find(DTemplate == 1 );
[ nonZPosY, nonZPosX ] = find( DTemplate == 1 );
Wsq = W.^2;
% wsqMat = {};
% for i = 1:size(DTemplate, 1)
%     fprintf('%d\n', i);
%     curX = find( DTemplate(i,:)~=0 );
%     for j = 1:length(curX)
%         wsqMat{curX(j), curX(j)}= ...
%             W( curX(j), :) .* W( curX(j), :);
%     end
%     for j = 1:length(curX)-1
%         for k = (j+1):length(curX)
%             wsqMat{curX(j), curX(k)} ...
%                 = W( curX(j), :) .* W( curX(k), :);
%         end
%     end
% end
nonZNum = 0;
for i = 1:size(DTemplate, 1)
    curX = find( DTemplate(i,:)~=0 );
    nonZNum = nonZNum + length( curX );
    for j = 1:length(curX)-1
        for k = (j+1):length(curX)
            nonZNum = nonZNum + 1;
        end
    end
end

%gennerate nonZYGrp, a group of need-to-update variables that have the same
%Y-value in the dictionary
%sorting by Y-vale in ascending order
%e.g. 1 4 6
%      3   7
%     2 5
%in this case, there are three groups {1,4,6} {3,7} {2,5}
[sNonZPosY, idx] = sort( nonZPosY );
PosGrp = zeros( length( unique( nonZPosY ) ), 2 );
cnt = 1;
for i = 1:size( sNonZPosY )
    target = sNonZPosY(i);
    for j = i+1:size( sNonZPosY )
        if sNonZPosY(j) ~= target
            break;
        end
    end
    PosGrp(cnt, :) = [i, j-1];
    cnt = cnt + 1;
end
nonZYGrp = [];
cnt = 1;
for i = 1:size( PosGrp, 1 )
    if PosGrp(i, 1) ~= PosGrp(i, 2)
        curGrp = idx(PosGrp(i, 1):PosGrp(i, 2));
        nonZYGrp{cnt} = curGrp;
        cnt = cnt + 1;
    end
end
%generating nonZGrpInfo, nonZLen, for adding up need-to-update variable
%nonZGrpInfo, a reverse map. Mapping from non-zero position to its
%dictionary element
%nonZLen, a vector [mLen, 1], with eacn entry records # non-zero variables
nonZGrpInfo = zeros( length( nonZPos ), 1 );
nonZLen = zeros( mLen, 1 );
curLoc = 1;
for i = 1:mLen
    len = length( find( DTemplate(:, i ) == 1 ) );
    nonZLen(i) = len;
    nonZGrpInfo(curLoc:(curLoc+len-1)) = i;
    curLoc = curLoc + len;
end
%generating constraints Hessian Matrix
% idx = 1:length(nonZPos);
% DConHes = sparse( idx, idx, 2, length(nonZPos), length(nonZPos) );
%assign the Hessian Function
Hes = @(curD, lambda)D_HessianFFunc(curD, lambda, W, Wsq, staticTerm, nonZPosY, nonZPosX, nonZPos, nonZLen, nonZYGrp, nonZNum, scaleFac);
% Hes = @(curD, lambda)D_HessianFFunc2(curD, lambda, W, staticTerm, wsqMat, nonZPosY, nonZPosX, nonZPos, nonZLen, nonZYGrp );
options.HessFcn = Hes;
lowerB = zeros( length( rD(nonZPos) ), 1 );
upperB = zeros( length( rD(nonZPos) ), 1 ); upperB(:) = Inf;

startD = full( rD(nonZPos) );
% objecFun = @(startD) D_termFunc( Y, startD, W, staticTerm, nonZPos, nonZPosY, nonZPosX, phi );
objecFun = @(startD) D_termFuncScale( Y, startD, W, staticTerm, nonZPos, nonZPosY, nonZPosX, phi, scaleFac );
% keyboard();
conFun = @(x)nlc_fn( x, nonZGrpInfo, nonZLen );
uD = fmincon( objecFun, startD, [], [], [], [], lowerB, upperB, conFun, options );
rD(nonZPos) = uD;

end
function [ val, grad ] = D_termFunc( inY, curD, curW, staticTerm, nonZPos, nonZPosY, nonZPosX, phi )
D = sparse( nonZPosY, nonZPosX, curD, size(inY, 1), size(curW, 1) );
preY = staticTerm + D * curW;
ePreY = exp( preY );
val = -inY.* preY + ePreY;
val = sum( sum( val ) );
val = val + phi* sum(curD);
grad = ( -inY + ePreY )*curW' + phi;
grad = grad( nonZPos );
end

function [ val, grad ] = D_termFuncScale( inY, curD, curW, staticTerm, nonZPos, nonZPosY, nonZPosX, phi, scaleFac)
D = sparse( nonZPosY, nonZPosX, curD, size(inY, 1), size(curW, 1) );
preY = staticTerm + D * curW;
ePreY = exp( preY );
val = -inY.* preY + ePreY;
val = sum( sum( val * scaleFac ) );
val = val + phi* sum(curD);
grad = ( ( -inY + ePreY ) * scaleFac )*curW' + phi;
grad = grad(nonZPos);
end

function [ c,ceq, GC, GCeq ]=nlc_fn( x, nonZGrpInfo, nonZLen )
%the number of constraints is the number of molecule species 
tmp = ( x.^2 );
sumLen = accumarray( nonZGrpInfo, tmp );
c =  sumLen - 1;
ceq = [];
GC = zeros( length( x ), length( nonZLen ) );
tmp = 2*x;
curLoc = 1;
for i = 1:length( nonZLen ) 
    GC(curLoc:(curLoc+nonZLen(i)-1), i) = tmp(curLoc:(curLoc+nonZLen(i)-1));
    curLoc = curLoc + nonZLen(i);
end
GC = sparse( GC );
GCeq= [];
end
function [ hMat ] = D_HessianFFunc( curD, lambda, curW, curWsq, staticTerm, nonZPosY, nonZPosX, nonZPos, nonZLen, nonZYGrp, nonZNum, scaleFac )
%setting the diagonal value of Hessian Matrix
D = sparse( nonZPosY, nonZPosX, curD, size( staticTerm, 1 ), size( curW, 1 ) );
preY = staticTerm + D * curW;
ePreY = exp( preY ) * scaleFac;
tmp = ( ePreY )*curWsq';
tmp = tmp(nonZPos);
HesLen = length( nonZPos );
hMat = sparse( [], [], [], HesLen, HesLen, nonZNum );
for i = 1:length( nonZPos )
    hMat(i, i)=tmp(i);
end
% setting the off-diagonal value of Hessian Matrix
% for i = 1:length( nonZYGrp )
%     curGrp = nonZYGrp{i};
%     grpY = nonZPosY(curGrp(1));
%     grpX = nonZPosX(curGrp );
%     for j = 1:(length( curGrp )-1)
%         tmpVal = ePreY(grpY, :).*curW(grpX(j), :);
%         tmpValAry = zeros( length( curGrp ) - j, 1 );
%         for z = j+1:length( curGrp )
%             tmpValAry(z-j) = tmpVal*curW( grpX(z), :)';
%         end
%         for z = j+1:length( curGrp )
%             hMat( curGrp(j), curGrp(z) ) = tmpValAry(z-j);
%             hMat( curGrp(z), curGrp(j) ) = tmpValAry(z-j);
%         end
%     end
% end
for i = 1:length( nonZYGrp )
    curGrp = nonZYGrp{i};
    grpY = nonZPosY(curGrp(1));
    grpX = nonZPosX(curGrp );
%     tmpVal = repmat( ePreY(grpY, :), length( curGrp ) - 1, 1 ).*curW( grpX(1:end-1), : );
    for j = 1:(length( curGrp )-1)
%         tmpValAry = zeros( length( curGrp ) - j, 1 );
%         curVal = tmpVal(j,:);
        curVal = ePreY(grpY, :).*curW(grpX(j),:);
        curGW = curW( grpX(j+1:length( curGrp )), : );
        for z = 1:(length( curGrp ) - j )
            tmpValAry(z) = curVal*curGW( z, :)';
        end
        for z = j+1:length( curGrp )
            hMat( curGrp(j), curGrp(z) ) = tmpValAry(z-j);
            hMat( curGrp(z), curGrp(j) ) = tmpValAry(z-j);
%             tmpVal = curVal*curW(grpX(z),:)';
%             hMat( curGrp(j), curGrp(z) ) = tmpVal;
%             hMat( curGrp(z), curGrp(j) ) = tmpVal;
        end
    end
end


curLoc = 1;
tmp = zeros(length(nonZPosY), 1);
for i = 1:length(nonZLen)
    tmp(curLoc:(curLoc+nonZLen(i)-1)) = lambda.ineqnonlin(i);
    curLoc = curLoc + nonZLen(i);
end
idx = 1:length(nonZPosY);
lambdaMatrix = sparse( idx, idx, tmp, length(nonZPosY), length(nonZPosY) );
lambdaMatrix = lambdaMatrix * 2;
hMat = sparse( hMat + lambdaMatrix );
end

%% Hessian Function with square Matrix of W precomputed
function [ hMat ] = D_HessianFFunc2( curD, lambda, curW, staticTerm, wsqMat, nonZPosY, nonZPosX, nonZPos, nonZLen, nonZYGrp )
%setting the diagonal value of Hessian Matrix
D = sparse( nonZPosY, nonZPosX, curD, size( staticTerm, 1 ), size( curW, 1 ) );
preY = staticTerm + D * curW;
ePreY = exp( preY );
HesLen = length( nonZPos );
hMat = sparse( HesLen, HesLen );

for i = 1:length( nonZYGrp )
    curGrp = nonZYGrp{i};
    grpY = nonZPosY(curGrp(1));
    grpX = nonZPosX(curGrp );
    for j = 1:(length( curGrp ))
        hMat( curGrp(j), curGrp(j) ) = ePreY(grpY, :)*wsqMat{grpX(j), grpX(j) }';
    end
    for j = 1:(length( curGrp )-1)
        for z = j+1:length( curGrp )
            tmp = ePreY(grpY, :)*wsqMat{grpX(j), grpX(z) }';
            hMat( curGrp(j), curGrp(z) ) = tmp;
            hMat( curGrp(z), curGrp(j) ) = tmp;
        end
    end
end

curLoc = 1;
tmp = zeros(length(nonZPosY), 1);
for i = 1:length(nonZLen)
    tmp(curLoc:(curLoc+nonZLen(i)-1)) = lambda.ineqnonlin(i);
    curLoc = curLoc + nonZLen(i);
end
idx = 1:length(nonZPosY);
lambdaMatrix = sparse( idx, idx, tmp, length(nonZPosY), length(nonZPosY) );
lambdaMatrix = lambdaMatrix * 2;
hMat = sparse( hMat + lambdaMatrix );
end