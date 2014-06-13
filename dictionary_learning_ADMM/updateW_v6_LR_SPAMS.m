function [ rW, rW0 ] = updateW_v6_LR_SPAMS( inY, outD, W_init, W0_init, rMatrix, DTemplate, lambda, theta)
%updateW_v6_LG_SPAMS using SPAMS package to solve the linear regression
%dictioanry learning problem, update W with ridge panelty on location and
%lasso term on W
param.lambda=lambda;
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
% and uses all the cores of the machine
param.mode=2; %
param.pos = true;
rW = W_init(:, :);
rW0 = W0_init(:);
[mLen, gLen] = size(rW);
Y = inY(:, :);
Y = log(Y); Y(Y==-inf)=0;
[sLen] = size(Y, 1);

preY = outD*rW + repmat( rW0', sLen, 1 );
preWa = [rW; rW0'];
threshold = 1e-6;
resY = Y - preY;
thetaRT = sqrt(theta);
rMatrix = rMatrix*thetaRT;
rLen = size( rMatrix, 1 );
resDiff = -rMatrix*rW';

opts.mFlag=0;       % treating it as compositive function 
opts.lFlag=0;       % Nemirovski's line search
opts.init=2;        % starting from a zero point
% termination criterion
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations

% normalization
opts.nFlag=0;       % without normalization

% regularization
opts.rFlag=1;    

%construct mD
xLocCel = cell( mLen, 1 );
yLocCel = cell( mLen, 1 );
dValCel = cell( mLen, 1 );
for i = 1:mLen
    fprintf('%d\n', i);
    spec =  DTemplate(:, i) ~= 0 ;
    curD = outD(spec, i);
    cSLen = length( curD );
    stPos = 1;
    nRow = cSLen*gLen;
    xLoc = zeros( nRow, 1);
    yLoc = zeros( nRow, 1);
    dVal = zeros( nRow, 1);
    % mD = zeros( sLen*gLen, gLen );
    for j = 1:gLen
        xLoc(stPos:(stPos+cSLen-1)) = j;
        yLoc(stPos:(stPos+cSLen-1)) = stPos:(stPos+cSLen-1);
        dVal(stPos:(stPos+cSLen-1)) = curD;
        stPos = stPos + cSLen;
    end
    xLocCel{i} = xLoc;
    yLocCel{i} = yLoc;
    dValCel{i} = dVal;
end
itNum = 1;
while 1
    fprintf('itNum = %d\n', itNum );
    for i = 1:mLen
%         fprintf('mLen = %d\n', i );
        curW = rW(i, :)';
        spec =  DTemplate(:, i) ~= 0 ;
        curD = outD(spec, i);
        [ mD ] = consMD( gLen, xLocCel{i}, yLocCel{i}, dValCel{i} );
        mD = [mD; rMatrix];
        resY(spec, :) = resY(spec, :) + ( curD*curW');
        cResY = resY(spec, :);
        resDiff(:,i) = resDiff(:,i) + (rMatrix*curW);
        cResDiff = resDiff(:, i);
        cResY = [cResY(:); cResDiff];
%         alpha = mexLasso( cResY, mD, param);
        [alpha, funVal2, ValueL2]= LeastR(mD, cResY, lambda, opts);
%         [alpha,FitInfo] = lasso(mD,cResY, 'Lambda', lambda );
        alpha = max( alpha, 0);
        rW(i,:) = alpha;
        resY(spec, :) = resY(spec, :) - (curD*alpha');
        resDiff(:,i) = resDiff(:,i) - (rMatrix*curW);
    end
    resY = resY + repmat(rW0', sLen, 1);
    rW0 = mean(resY, 1)';
    resY = resY - repmat(rW0', sLen, 1);
    toc
    rWa = [rW; rW0'];
    fprintf( '%g %g\n', max( abs( rWa(:) - preWa(:) ) ), threshold);
    if  max( abs( rWa(:) - preWa(:) ) ) <= threshold
        break;
    end
    itNum = itNum + 1;
    preWa = rWa;
end
rW = rWa(1:(end-1), :);
rW0 = rWa(end, :);
end

function [ mD ] = consMD( gLen, xLoc, yLoc, dVal )
% sLen = length( dVec );
% stPos = 1;
% nRow = sLen*gLen;
% xLoc = zeros( nRow, 1);
% yLoc = zeros( nRow, 1);
% dVal = zeros( nRow, 1);
% % mD = zeros( sLen*gLen, gLen );
% for i = 1:gLen
%     xLoc(stPos:(stPos+sLen-1)) = i;
%     yLoc(stPos:(stPos+sLen-1)) = stPos:(stPos+sLen-1);
%     dVal(stPos:(stPos+sLen-1)) = dVec;
% %     mD(stPos:(stPos+sLen-1),i) = dVec;
%     stPos = stPos + sLen;
% end
nRow = length( xLoc );
mD = sparse( yLoc, xLoc, dVal,  nRow, gLen );
end

function [val] = ll( inY, inD, inW, inW0, lambda, theta)
    val = sum((inY - inD*inW - inW0).^2);
    val = val + lambda*inW
end