function [ rD ] = updateD_v3( inY, outW, outW0, D_init, DTemplate, aMatrix, phi  )
%updateD_v3 update one dictionary element at a time using fminunc without ADMM
[sLen, mLen] = size(D_init);

options = optimoptions('fminunc');
options.MaxFunEvals= 1e6;
options.MaxIter = 100;
options.TolFun=1e-6;
% options.TolCon=1e-8;
% options.Algorithm='interior-point';
options.Algorithm='trust-region';
options.Display='none';
% options.AlwaysHonorConstraints = 'none';
% options.UseParallel='always';
options.GradObj='on';
% options.GradConstr = 'on';
Y = inY(:, :);
Y = Y(:, aMatrix == 1);
W = outW(:, :);
W = W(:,  aMatrix == 1  );
outW0 = outW0( aMatrix == 1 );
preY = D_init*W + repmat( outW0(:)', sLen, 1 );
rD = D_init;
prevD = D_init;
Wsq = W.^2;
threshold = 1e-4;
options.Hessian='on';
WIdx = cell( mLen, 1 );
WsqIdx = cell( mLen, 1 );
for i = 1:mLen
    WIdx{i} = W(i, :);
    WsqIdx{i} = Wsq(i, :);
end

while 1
    tic
    for i = 1:mLen
%         fprintf('%d\n', i);
%         curW = W(i, :);
%         curWsq = Wsq(i,:);
        nonZPos = find(DTemplate(:,i)~=0);
        resPreY = preY(nonZPos, :) - rD(nonZPos, i)*WIdx{i};
%         lowerB = zeros( length(nonZPos), 1 );
%         upperB = zeros( length(nonZPos), 1 ); upperB(:) = Inf;
        test = @(curD) D_termFunc( Y(nonZPos, :), curD, WIdx{i}, resPreY, phi, WsqIdx{i} );
            %the last two lines only valid for interior-point algo.
%             options.Hessian = 'user-supplied';
%             options.HessFcn = @(curD, lambda)D_HessianFFunc(curD, lambda, curW, curWsq, resPreY);
%         [uD,~] = fmincon(test,rD(nonZPos, i),[],[],[],[],lowerB,upperB, @nlc_fn, options);

%         options.Hessian='on';
        uD = fminunc(test,rD(nonZPos, i),options);
        uD = max(uD, 0);
% % %         
        uD = uD/(max(1, norm( uD, 2 ) ) );
        rD(nonZPos, i) = uD;
        preY(nonZPos, :) = resPreY + uD*WIdx{i};
        %         fprintf( '%g\n', testLP1 - testLP2 );
%         if testLP1 - testLP2 < 0 && abs(testLP1-testLP2) > 1e-6
%             keyboard();
%         end
    end
    toc
    fprintf( '%g %g\n', max( abs( rD(:) - prevD(:) ) ), threshold);
    if  max( abs( rD(:) - prevD(:) ) ) <= threshold
        break;
    end
    prevD = rD;
end

end
function [ val, grad, hMat ] = D_termFunc( inY, curD, curW, resPreY, phi, curWsq )
%inY: [s(non-zero-entry) a*b]
%curD: [s(non-zero-entry) 1]
%curW: [1 a*b]
%resPreY: [s(non-zero-entry) a*b]
preY = resPreY + curD*curW;
ePreY = exp( preY );
val = -inY.*preY + ePreY;
val = sum(val(:)) + phi * sum( curD );
grad = (-inY+ePreY)*curW' + phi;
tmp = (ePreY)*curWsq';
sLen = size(curD, 1);
idx = 1:sLen;
hMat = sparse( idx, idx, tmp, sLen, sLen );
end
function [c,ceq, gradc, gradceq]=nlc_fn(x)
c = sum(x.^2)- 1;
ceq = [];
gradc = 2*x;
gradceq = [];
end
function [hMat] = D_HessianFFunc(curD, lambda, curW, curWsq, resPreY)
preY = resPreY + curD*curW;
ePreY = exp( preY );
tmp = (ePreY)*curWsq';
sLen = size(curD, 1);
idx = 1:sLen;
hMat = sparse( idx, idx, tmp, sLen, sLen );
lambdaMatrix = sparse( idx, idx, 2*lambda.ineqnonlin, sLen, sLen );
hMat = sparse( hMat + lambdaMatrix );
end
