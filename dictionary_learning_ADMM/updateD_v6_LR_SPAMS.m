function [ rD ] = updateD_v6_LR_SPAMS( inY, outW, outW0, D_init, DTemplate, aMatrix, phi)
%updateD_v6_LG_SPAMS using SPAMS package to solve the linear regression
%dictioanry learning problem
%   Detailed explanation goes here
param.lambda=phi; % not more than 20 non-zeros coefficients
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
% and uses all the cores of the machine
param.mode=2; %
param.pos = true;
[sLen, mLen] = size(D_init);
Y = inY(:, :);
Y = Y(:, aMatrix == 1);
Y = log(Y); Y(Y==-inf)=0;
Y = Y';
W = outW(:, :);
W = W(:,  aMatrix == 1  );
gLen = size( W, 2 );
outW0 = outW0(aMatrix==1);
preY = D_init*W + repmat( outW0(:)', sLen, 1 );
preY = preY';
rD = D_init;
prevD = D_init;
threshold = 1e-6;
resY = Y - preY;

WIdx = cell( mLen, 1 );
resYIdx = cell( mLen, 1 );
for i = 1:mLen
    WIdx{i} = W(i,:);
    nonZPos = DTemplate(:,i)~=0;
end
while 1
    tic
    for i = 1:mLen
%         curW = W(i, :);
        nonZPos = find(DTemplate(:,i)~=0);
        [ mW ] = consMW( WIdx{i}, length( nonZPos ) );
%         mW = full( mW );
        curD = rD( nonZPos, i );
        cResY = resY(:,nonZPos);
        cResY = cResY + (curD*WIdx{i})';
%         curY = resY(:,nonZPos);
        alpha=mexLasso( cResY(:),mW,param);
        alpha = alpha / ( max( norm( alpha, 2 ), 1 ) );
        rD(nonZPos, i) = alpha;
        cResY = cResY - (alpha*WIdx{i})';
        resY(:, nonZPos ) = cResY;
    end
    toc
    fprintf( '%g %g\n', max( abs( rD(:) - prevD(:) ) ), threshold);
    if  max( abs( rD(:) - prevD(:) ) ) <= threshold
        break;
    end
    prevD = rD;
end
end

function [ mW ] = consMW( wVec, sLen )
%mW initialize only one times
gLen = length( wVec );
% nrow = gLen*sLen;
% xLoc = zeros( nrow, 1 );
% yLoc = zeros( nrow, 1 );
% wVal = zeros( nrow, 1 );
sPos = 1;
% mW = gmW( 1:(gLen*sLen), 1:sLen );
% mW(:) = 0;
mW = zeros( gLen*sLen, sLen );
for i = 1:sLen
%     xLoc(sPos:(sPos+gLen-1)) = i;
%     yLoc(sPos:(sPos+gLen-1)) = sPos:(sPos+gLen-1);
%     wVal(sPos:(sPos+gLen-1)) = wVec;
    mW(sPos:(sPos+gLen-1),i) = wVec;
    sPos = sPos + gLen;
end
% mW = sparse( yLoc, xLoc, wVal, gLen*sLen, sLen );
end
