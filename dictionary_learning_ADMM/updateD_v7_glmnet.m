function [ uD ] = updateD_v7_glmnet( inY, outW, outW0, D_init, DTemplate, aMatrix, phi, itNum )
%updateD_v7_glmnet update the whole dictionary using fmincon without ADMM
%aMatrix, a indicator matrix, with 1 means using in trainning and 0 means
%using in testing
%phi, the sparsity weighting parameter
%Error occurred when updating with user-defined lambda sequence,
%deprecated.
%% setting glmnet parameters
options = glmnetSet;
options.nlambda=1;
options.lambda=[phi];
options.cl = [0; inf];
options.alpha = 1;
options.maxit = 1e5;
options.intr = true;
options.standardize = false;
%% setting input data
[sLen, mLen] = size(D_init);
Y = inY(:, :);
Y = Y(:, aMatrix == 1);
W = outW(:, :);
W = W(:,  aMatrix == 1 );
ins = zeros( mLen, 1 );
for i = 1:mLen
    ins(i) = max( W(i, :) );
end
mIdx = find( ins>1e-3 );
rMLen = length(mIdx);
W = W(mIdx, :);
nLen = size( W, 2 );
outW0 = outW0(aMatrix==1);
rD = D_init(:, mIdx);
uD = D_init;
prevD = rD;
preY = rD*W + repmat( outW0(:)', sLen, 1 );
threshold = 1e-4;

%sLen collections
% sLenCol = [];
% for i = 1:mLen
%     nonZPos = find(DTemplate(:,i)~=0);
%     sLenCol = [ sLenCol length( nonZPos ) ];
% end
% sLenCol = unique( sLenCol );
%construct W collection
WCol = cell( rMLen, 1 );
for i = 1:rMLen
    curW = W(i,:);
    nonZPos = find(DTemplate(:,i)~=0);
    curSLen = length( nonZPos );
    curWVec = repmat( curW, curSLen, 1 );
    sVec = curWVec(:);
    iSize = curSLen*nLen;
    jSize = curSLen;
    jCoord = (1:curSLen)';
    jCoord = repmat( jCoord, nLen, 1 );
    jCoord = jCoord(:);
    iCoord = 1:iSize;
    WCol{i} = sparse(iCoord, jCoord, sVec, iSize, jSize, iSize );
end
%construct w index
WIdx = cell( rMLen, 1 );
for i = 1:rMLen
    WIdx{i} = W(i,:);
end

while 1
    tic
    for i = 1:rMLen
%         fprintf('%d\n', i );
%         curW = W(i, :);
        nonZPos = find(DTemplate(:,i)~=0);
        curD = rD(nonZPos, i);
        resPreY = preY(nonZPos, :) - curD*WIdx{i};
        tmpBeta = zeros(length(nonZPos), 1);
        for j = 1:length(nonZPos)
            curW = W(i,:)';
            curY = Y(nonZPos(j), :);
            curY = curY(:);
            tmp = resPreY(j, :); 
            options.offset = tmp(:); 
            fit = glmnet( curW, curY, 'poisson', options);
            %max(curY'-tmp)
             tmpBeta(j) = fit.beta(end);
%             [B] = lassoglm(curW, curY, 'poisson', 'Lambda', phi, 'Offset', tmp, 'Standardize', false);
            tmpBeta(j) = fit.beta(1);
        end
%         B = max( B, 0 );
        tmpBeta = tmpBeta / max( norm(tmpBeta), 1 );
        rD(nonZPos, i) = tmpBeta;
        preY(nonZPos, :) = resPreY + rD(nonZPos, i)*WIdx{i};
%         curW = WCol{i};
%         curD = rD(nonZPos, i);
%         resPreY = preY(nonZPos, :) - curD*WIdx{i};
%         curY = inY(nonZPos,:);
%         curY = curY(:);
%         options.offset = resPreY(:);
%         fit = glmnet( curW, curY, 'poisson', options);
% %             [B,FitInfo] = lassoglm( curW, curY, 'poisson', 'Standardize', false, 'Lambda', phi, 'Offset', cResY(:) ); 
%         B = fit.beta(:,1);
%         B = B / max( norm(B), 1 );
%         rD(nonZPos, i) = B;
%         preY(nonZPos, :) = resPreY + rD(nonZPos, i)*WIdx{i};
    end
    toc
    fprintf( '%g %g\n', max( abs( rD(:) - prevD(:) ) ), threshold);
    if  max( abs( rD(:) - prevD(:) ) ) <= threshold
        break;
    end
    prevD = rD;
end
uD(:, mIdx) = rD;
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