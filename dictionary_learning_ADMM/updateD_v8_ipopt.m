function [ finalD ] = updateD_v8_ipopt( inY, outW, outW0, D_init, DTemplate, aMatrix, HesOpt, phi, scaleFac, itNum )
%updateD_v8 update the whole dictionary using fmincon without ADMM
%aMatrix, a indicator matrix, with 1 means using in trainning and 0 means
%using in testing
%HesOpt, a flag of using Hessian or not, 1 means using. deprecated, not
%used in this implementation
%phi, the sparsity weighting parameter
%scaleFac, scalar, to make phi more controlable

%% initialize variables
[sLen, mLen] = size(D_init);
aIdx = aMatrix(:)==1;
Y = sparse(inY(:, aIdx));
W = outW(:, aIdx);
staticTerm = repmat( outW0(aIdx)', sLen, 1 );
ins = zeros( mLen, 1 );
for i = 1:mLen
    ins(i) = max( W(i, :) );
end
rMIdx = find( ins > 1e-3 );
rMLen = length( rMIdx );
W = W( rMIdx, : );
rD = full( D_init( :, rMIdx ) );
DTemplate = DTemplate(:, rMIdx );
nonZPos = find( DTemplate == 1 );
[ nonZPosY, nonZPosX ] = find( DTemplate == 1 );

%% option of ipopt section
% The starting point.
startD = rD(nonZPos);       
% The constraint functions are bounded from below by zero.
options.lb = zeros(length(nonZPos),1);
options.ub = ones(length(nonZPos),1);
options.cl = -Inf(rMLen,1);
options.cu = ones(rMLen,1);
% Set up the auxiliary data.

% Set the IPOPT options.
options.ipopt.print_level = 5;
% options.ipopt.jac_d_constant   = 'yes';
% options.ipopt.hessian_constant = 'yes';
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         = itNum;
options.ipopt.tol              = 1e-8;
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.limited_memory_max_history = 1500; 
options.ipopt.honor_original_bounds = 'yes';
options.ipopt.bound_relax_factor = 0;

% The callback functions.
funcs.objective         = @objective;
funcs.constraints  = @constraints;
funcs.gradient          = @gradient;
funcs.jacobian          = @jacobian;
funcs.jacobianstructure = @jacobianstructure;
%   
% funcs.hessian           = @hessian;
% funcs.hessianstructure  = @hessianstructure;

%% generating metadata for nonZeroPosition
if HesOpt == 1
    fprintf( 'No computing Hessian of D, this options is currently under construction...\n' );
end
%generate nonZPos (one-dimension indicator), map nonzero element -> element
%in the dictionary. nonZPosY (two-dimension indicator, y-axis), map nonzero
%element -> element's Y position in the dictionary
%nonZPosX (two-dimension indicator, x-axis) map nonzero
%element -> element's X position in the dictionary
%nonZPos = find(DTemplate == 1 );
% Wsq = W.^2;
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
% nonZNum = 0;
% for i = 1:size(DTemplate, 1)
%     curX = find( DTemplate(i,:)~=0 );
%     nonZNum = nonZNum + length( curX );
%     for j = 1:length(curX)-1
%         for k = (j+1):length(curX)
%             nonZNum = nonZNum + 1;
%         end
%     end
% end

%gennerate nonZYGrp, a mapping from y-axis to # variables
%e.g. nonZYGrp{1} = {2, 5} means in y-axis = 1, there are two variables,
%variables 2 and 5
% [sNonZPosY, idx] = sort( nonZPosY );
% PosGrp = zeros( length( unique( nonZPosY ) ), 2 );
% cnt = 1;
% for i = 1:size( sNonZPosY )
%     target = sNonZPosY(i);
%     for j = i+1:size( sNonZPosY )
%         if sNonZPosY(j) ~= target
%             break;
%         end
%     end
%     PosGrp(cnt, :) = [i, j-1];
%     cnt = cnt + 1;
% end
% nonZYGrp = [];
% cnt = 1;
% for i = 1:size( PosGrp, 1 )
%     if PosGrp(i, 1) ~= PosGrp(i, 2)
%         curGrp = idx(PosGrp(i, 1):PosGrp(i, 2));
%         nonZYGrp{cnt} = curGrp;
%         cnt = cnt + 1;
%     end
% end
%generating nonZGrpInfo, nonZLen, for adding up need-to-update variable
%nonZGrpInfo, a reverse map. Mapping from non-zero position to its
%dictionary element
%nonZLen, a vector [mLen, 1], with eacn entry records # non-zero variables
nonZGrpInfo = zeros( length( nonZPos ), 1 );
nonZLen = zeros( rMLen, 1 );
curLoc = 1;
for i = 1:rMLen
    len = length( find( DTemplate(:, i ) == 1 ) );
    nonZLen(i) = len;
    nonZGrpInfo(curLoc:(curLoc+len-1)) = i;
    curLoc = curLoc + len;
end
Wsq = W.^2;  
WT = W';
gRes = zeros( size( rD ) );
gRes = gRes( nonZPos );
gPreY = zeros( size(Y) );
gEPreY = zeros( size(Y) );
jacobStruct = zeros(rMLen, length( nonZPos ) );
for i = 1:rMLen
    target =  nonZPosX==i ;
    jacobStruct(i, target) = 1;
end
jacobStruct = sparse( jacobStruct );

nonZIdx = cell( rMLen, 1 );
for i = 1:rMLen
    nonZIdx{i} = find( DTemplate(:, i) > 0 );
end
YT=Y';
sTermT = staticTerm';
% HesLen = length( nonZPos );
% HessianStruct = zeros( HesLen, HesLen );
% for i = 1:HesLen
%     HessianStruct(i, i)=1;
% end
% for i = 1:length( nonZYGrp )
%     curGrp = nonZYGrp{i};
%     for j = 1:length(curGrp)-1
%         for z = j+1:length(curGrp)
%             HessianStruct(curGrp(j),  curGrp(z)) = 1;
%         end
%     end
% end

% options.auxdata = { Y, W, staticTerm, nonZPos, nonZPosY, nonZPosX, phi, scaleFac, nonZGrpInfo, nonZLen, jacobStruct, nonZYGrp, Wsq, HessianStruct };
options.auxdata = { Y, W, staticTerm, nonZPos, nonZPosY, nonZPosX, phi, scaleFac, nonZGrpInfo, nonZLen, jacobStruct, WT, gRes, nonZIdx, YT, sTermT, gPreY, gEPreY};
[x, ~] = ipopt_auxdata(startD,funcs,options);
rD(nonZPos) = x;
% for i = 1:rMLen
%     rD(:, i) = rD(:, i) / max( 1, norm( rD(:, i) ) ); 
% end
finalD = sparse( D_init );
finalD( :, rMIdx ) = rD;
end

function f = objective (x, auxdata)
%x is current D
[Y, W, ~, nonZPosY, nonZPosX, phi, scaleFac , WT, YT, sTermT] = deal(auxdata{[1 2 3 5 6 7 8 12 15 16]});
% D = sparse( nonZPosY, nonZPosX, x, size(Y, 1), size(W, 1) );
% gPreY = staticTerm + D * W;
% gEPreY = exp( gPreY );
DT = sparse( nonZPosX, nonZPosY, x, size(W, 1), size(Y, 1) );
preY = sTermT + WT * DT;
ePreY = exp( preY );
% f = -Y.* preY + ePreY;
f = sum( sum( -YT.* preY + ePreY ) ) * scaleFac;
f = f + phi* sum(x);
end

function gRes = gradient (x, auxdata)
[Y, W, ~, ~, nonZPosY, nonZPosX, phi, scaleFac, nonZLen, WT, gRes, nonZIdx, YT, sTermT] = deal(auxdata{[1:8 10 12 13 14 15 16]});
DT = sparse( nonZPosX, nonZPosY, x, size(W, 1), size(Y, 1) );
preY = sTermT + WT * DT;
ePreY = exp( preY );
% gRes = (-Y+ePreY)*WT;
% gRes = gRes(nonZPos);
curLoc = 1;
for i = 1:length(nonZIdx)
    idx = nonZIdx{i};
    gRes(curLoc:(curLoc+nonZLen(i)-1)) = (-YT(:,idx) + ePreY(:, idx))'*WT(:,i);
    curLoc = curLoc + nonZLen(i);
end
gRes = gRes*scaleFac + phi;
end

function c = constraints(x, auxdata)
[nonZGrpInfo]=auxdata{9};
tmp = ( x.^2 );
c = accumarray( nonZGrpInfo, tmp );
end

function J = jacobian (x, auxdata)  
nonZLen = auxdata{10};
J = zeros( length( nonZLen ), length( x ) );
tmp = 2*x;
curLoc = 1;
for i = 1:length( nonZLen ) 
    J(i, curLoc:(curLoc+nonZLen(i)-1)) = tmp(curLoc:(curLoc+nonZLen(i)-1));
    curLoc = curLoc + nonZLen(i);
end
J=sparse(J);
end

function J = jacobianstructure (auxdata)
[J] = deal(auxdata{11});
end

function H = hessianstructure (auxdata)
[H] = deal(auxdata{[end]});
end
  
function  H = hessian(x, sigma, lambda, auxdata)  
%setting the diagonal value of Hessian Matrix
[ W, staticTerm, nonZPos, nonZPosY, nonZPosX, scaleFac, nonZLen, nonZYGrp, Wsq, HessianStruct] = deal(auxdata{[2 3 4 5 6 8 10 12 13 14]});
D = sparse( nonZPosY, nonZPosX, x, size( staticTerm, 1 ), size( W, 1 ) );
preY = staticTerm + D * W;
ePreY = exp( preY ) * scaleFac;
tmp = ( ePreY )*Wsq';
tmp = tmp(nonZPos);
HesLen = length( nonZPos );
% H = sparse( [], [], [], HesLen, HesLen, length(find(HessianStruct==1)) );
H = zeros( HesLen, HesLen );
for i = 1:HesLen
    H(i, i)=tmp(i);
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
        curVal = ePreY(grpY, :).*W(grpX(j),:);
        curGW = W( grpX(j+1:length( curGrp )), : );
        for z = 1:(length( curGrp ) - j )
            tmpValAry(z) = curVal*curGW( z, :)';
        end
        for z = j+1:length( curGrp )
            H( curGrp(j), curGrp(z) ) = tmpValAry(z-j);
            H( curGrp(z), curGrp(j) ) = tmpValAry(z-j);
%             tmpVal = curVal*curW(grpX(z),:)';
%             hMat( curGrp(j), curGrp(z) ) = tmpVal;
%             hMat( curGrp(z), curGrp(j) ) = tmpVal;
        end
    end
end


curLoc = 1;
tmp = zeros(length(nonZPosY), 1);
for i = 1:length(nonZLen)
    tmp(curLoc:(curLoc+nonZLen(i)-1)) = lambda(i);
    curLoc = curLoc + nonZLen(i);
end
idx = 1:length(nonZPosY);
lambdaMatrix = sparse( idx, idx, tmp, length(nonZPosY), length(nonZPosY) );
lambdaMatrix = lambdaMatrix * 2;
H = sparse( tril( sigma*H + lambdaMatrix ) );
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