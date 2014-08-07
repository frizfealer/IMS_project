function [ uD ] = updateD_v3( inY, outW, outW0, D_init, DTemplate, aMatrix, phi  )
%updateD_v3 update one dictionary element at a time using fminunc or (minFunc) without ADMM
% aMatrix: the matrix with the same size of samples, indicating which
% samples to be used in fitting D (value = 1) and which are not (value = 0)
%deprecated, not applicable to L1-penalty
[sLen, mLen] = size(D_init);
[nLen] = size(outW, 2);

%% options for fminunc
options = optimoptions('fminunc');
options.MaxFunEvals= 1e6;
options.MaxIter = 100;
options.TolFun=1e-6;
options.Algorithm='trust-region';
options.Display='none';
options.GradObj='on';
options.Hessian='on';

%% options for minFunc
options=[];
options.display = 'none';
options.method='newton';

%% prepared variables
Y = inY(:, :);
Y = Y(:, aMatrix == 1);
W = outW(:, :);
W = W(:,  aMatrix == 1  );
ins = zeros( mLen, 1 );
for i = 1:mLen
    ins(i) = max( W(i, :) );
end
mIdx = find( ins>1e-3 );
rMLen = length(mIdx);
W = W(mIdx, :);
outW0 = outW0( aMatrix == 1 );
uD = D_init;
rD = D_init(:, mIdx);
preY = rD*W + repmat( outW0(:)', sLen, 1 );
prevD = rD;
Wsq = W.^2;

TolX = 1e-3;
TolFun = 1e-4;
%% construct W index and Y index
WIdx = cell( rMLen, 1 );
WsqIdx = cell( rMLen, 1 );
YIdx = cell( rMLen, 1 );
for i = 1:rMLen
    WIdx{i} = W(i, :);
    WsqIdx{i} = Wsq(i, :);
    nonZPos = DTemplate(:,i)~=0;
    YIdx{i} = Y(nonZPos, :);
end

prevFVal = D_termEval( Y, preY, phi, rD );
cnt = 1;
while 1
    tic
    for i = 1:rMLen
%         fprintf('%d\n', i);
%         curW = W(i, :);
%         curWsq = Wsq(i,:);
        nonZPos = find(DTemplate(:,i)~=0);
        tmpD = rD(nonZPos, i);
        resPreY = preY(nonZPos, :) - tmpD*WIdx{i};
%         lowerB = zeros( length(nonZPos), 1 );
%         upperB = zeros( length(nonZPos), 1 ); upperB(:) = Inf;
        test = @(curD) D_termFunc( YIdx{i}, curD, WIdx{i}, resPreY, phi, WsqIdx{i} );
            %the last two lines only valid for interior-point algo.
%             options.Hessian = 'user-supplied';
%             options.HessFcn = @(curD, lambda)D_HessianFFunc(curD, lambda, curW, curWsq, resPreY);
%         [uD,~] = fmincon(test,rD(nonZPos, i),[],[],[],[],lowerB,upperB, @nlc_fn, options);

%         options.Hessian='on';
%         uD = fminunc(test,rD(nonZPos, i),options);

        resD = minFunc(test,tmpD,options);
        resD = max(resD, 0);
        resD = resD/(max(1, norm( resD, 2 ) ) );
        rD(nonZPos, i) = resD;
        preY(nonZPos, :) = resPreY + resD*WIdx{i};
        %         fprintf( '%g\n', testLP1 - testLP2 );
%         if testLP1 - testLP2 < 0 && abs(testLP1-testLP2) > 1e-6
%             keyboard();
%         end
    end
    toc
    FVal = D_termEval( Y, preY, phi, rD );
    fprintf( 'it %d: TolX: %g %g, TolFun: %g %g, norm of step: %g\n',cnt, max( abs( rD(:) - prevD(:) ) ), TolX, FVal-prevFVal, TolFun, norm( rD(:)-prevD(:) ));
    if  max( abs( rD(:) - prevD(:) ) ) <= TolX || abs(FVal-prevFVal) <= TolFun
        break;
    end
    prevD = rD;
    prevFVal = FVal;
    cnt = cnt+1;
end
uD(:, mIdx) = rD;
end
function [ val, grad, hMat] = D_termFunc( inY, curD, curW, resPreY, phi, curWsq )
%inY: [s(non-zero-entry) a*b]
%curD: [s(non-zero-entry) 1]
%curW: [1 a*b]
%resPreY: [s(non-zero-entry) a*b]
curD = max( 0, curD );
preY = resPreY + curD*curW;
ePreY = exp( preY );
val = -inY.*preY + ePreY;
val = sum(val(:)) + phi * sum( abs(curD) );
% tmp = ones( size( curD ) )*phi;
% tmp(curD<0)=tmp(curD<0)*-1;
grad = (-inY+ePreY)*curW' + phi;
tmp = (ePreY)*curWsq';
sLen = size(curD, 1);
idx = 1:sLen;
hMat = sparse( idx, idx, tmp, sLen, sLen );
end
% function [c,ceq, gradc, gradceq]=nlc_fn(x)
% c = sum(x.^2)- 1;
% ceq = [];
% gradc = 2*x;
% gradceq = [];
% end
% function [hMat] = D_HessianFFunc(curD, lambda, curW, curWsq, resPreY)
% preY = resPreY + curD*curW;
% ePreY = exp( preY );
% tmp = (ePreY)*curWsq';
% sLen = size(curD, 1);
% idx = 1:sLen;
% hMat = sparse( idx, idx, tmp, sLen, sLen );
% lambdaMatrix = sparse( idx, idx, 2*lambda.ineqnonlin, sLen, sLen );
% hMat = sparse( hMat + lambdaMatrix );
% end


function [ val] = D_termEval( inY, preY, phi, curD )
%inY: [s(non-zero-entry) a*b]
%curD: [s(non-zero-entry) 1]
%curW: [1 a*b]
%resPreY: [s(non-zero-entry) a*b]
ePreY = exp( preY );
val = -inY.*preY + ePreY;
val = sum(val(:)) + phi * sum( abs(curD(:)) );
end