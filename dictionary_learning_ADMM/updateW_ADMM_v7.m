function [A] = updateW_ADMM_v6( InputData, D, A, verbose, varargin )
%--------------------------------------------------------------------------
%updateW_ADMM_v6: update the AdmmModel based on inputData and D
%--------------------------------------------------------------------------
% DESCRIPTION:
%
% INPUT ARGUMENTS: 
%   InputData, the InputData data structure.
%   D, the dictionary
%   A, the AdmmModel data structure, can be initialize by initWADMM
%   verbose, whether to have verbose output
%   varargin, additional parameters, under construction
% OUTPUT ARGUMENTS: 
%   A, the updated AdmmModel data structure


%% correcting RHO_INIT when the W is too deviate from the data.
if max(A.W(:)) == 0 
    %in the initialzation of W, if hyperparameters == 0
    %we set RHO_INIT to be at most 1e-8.
    if A.RHO_INIT > 1e-8 && ...
            ( A.lambda == 0 && A.theta == 0 )
        A.RHO_INIT = 1e-8;
    end
    %in the initialzation of W, if one of the hyperparameter ~= 0
    %we set RHO_INIT to be a least 1e-6.
    if A.RHO_INIT < 1e-6 && ...
            ( A.lambda ~= 0 || A.theta ~= 0 )
        A.RHO_INIT = 1e-6;
    end
end

%% pour out data
Y = InputData.dataCube;
BlkDS = InputData.BlkDS;
[sLen, ~, ~] = size( Y );
aMatrix = InputData.aMatrix;
Rblk = InputData.Rblk;
scaleFactor = InputData.scaleFactor;
logFY = InputData.logFY;

%% variable setting
newWInfo = [];
if ~isempty(varargin)
    if length(varargin) >=1
        newWInfo = varargin{1};
    end
    if length(varargin) >= 2
        kappa = varargin{2};
    end
end

if verbose == 1
    fprintf( 'iteration number for W_update_ADMM = %d\n', A.itNum );
    temp = (aMatrix==0)&(BlkDS.indMap==1);
    fprintf( 'held-out sample number = %d, training sample number = %d\n', length( find( temp == 1 ) ), length( find( aMatrix == 1 ) ) );
    fprintf( 'lambda = %d, theta = %d\n', A.lambda, A.theta );
end
%% modifying D, if a dictionary element is too small, remove it.
ins = max( D, [], 1 );
usedMolIdx = find( ins > A.D_LOWER_BOUND );
D = D(:, usedMolIdx);
mLen = length( usedMolIdx );

%% additional fucnction for held-out validation, currently under construction
if ~isempty( newWInfo )
    for i = 1:length(newWInfo)
        if BlkDS.indMap(i) == 1
            newWInfo{i} = intersect( newWInfo{i}, usedMolIdx );
        end
    end
end

%% setting u0, u1, u2 and
u0 = sparse( zeros( size( A.z0(:, :) ) ) );
u1 = sparse( zeros( size( A.z1(usedMolIdx, :) ) ) );
u2 = cell( BlkDS.blkNum, 1 );
for i = 1:BlkDS.blkNum
    u2{i} = sparse( zeros( size( A.z2{i}(:, usedMolIdx) ) ) );
end

%% initialize others
LPAryADMM = zeros( A.itNum+1, BlkDS.blkNum );
resRecAry = zeros( A.itNum, 4, BlkDS.blkNum );
curRhoAry = zeros( A.itNum, BlkDS.blkNum );
curRhoAry(1, :) = A.RHO_INIT;


moutD = [D ones(sLen, 1)];
tmpEye = eye(mLen+1, mLen+1);
tmpEye(end, end) = 0;
aPartofCforUW = sparse( [moutD; tmpEye] );
stuckFlag = 0;
oITNum = A.itNum;
%% main optimization process
for j = 1:BlkDS.blkNum %for a block
    LPAryADMM(1, j) = LP_DL_Poiss( A.LINK_FUNC, aMatrix, Y, ...
    A.W(usedMolIdx, :), A.W0, D, A.lambda, 0, A.theta, scaleFactor, logFY, 0 );
    %initialize variables in a block
    loc = BlkDS.B2GMap{j};
    curY = Y(:, loc); curW = A.W(usedMolIdx, loc); curW0 = A.W0(loc);
    curZ0 = A.z0(:, loc); curU0 = u0(:, loc);
    curZ1 = A.z1(usedMolIdx, loc); curU1 = u1(:, loc);
    curZ2 = A.z2{j}(:, usedMolIdx); curU2 = u2{j};
    cRblk = Rblk{j};
    CforUW = replicateC( aPartofCforUW, length(loc) );
    RforUW = genSparseGroupingMatrix2( cRblk, mLen, 1 );
    CforUW = rowCombSparseMatrix( CforUW, RforUW );
    %%%%%additional fucnction for held-out validation, currently under construction
    nrPart1 = size(CforUW, 1);
    if ~isempty(newWInfo)
        [CforUW, ~] =  modifyCWithSparsity( CforUW, newWInfo, loc, mLen+1, nrPart1 );
    end
    %%%%%
    if verbose == 1
        fprintf( 'LP %g\n', LPAryADMM(1, j) );
    end
    itNumADMM = 1;
    while itNumADMM <= A.itNum
        preW = curW;
        curRho = curRhoAry(itNumADMM, j);
        %% update W
        if verbose == 1
            fprintf( 'updating W... ' );
        end
        [ w_w0 ] = update_w_w0_v4( curZ0, curU0, curZ1, curU1, curZ2, curU2, CforUW, curRho, [curW; curW0(:)'], 1 );
        tmp = reshape(w_w0, mLen+1, length(loc));
        curW = tmp(1:(end-1), :);
        curW0 = tmp(end,:);
        %% update Z0
        if verbose == 1
            fprintf( 'updating z0... ' );
        end
        locAlpha0 = intersect( loc, find( aMatrix(:) == 0 ) );
        locAlpha0 = find( ismember( loc, locAlpha0 ) );
        locAlpha1 = intersect( loc, find( aMatrix(:) == 1 ) );
        locAlpha1 = find( ismember( loc, locAlpha1 ) );
        %update for held-out data
        preZ0 = curZ0;
        [ curZ0(:, locAlpha0) ] = updatez0_v4( 0, A.LINK_FUNC, curY(:, locAlpha0), curZ0(:, locAlpha0), D, ...
            curW(:, locAlpha0), curW0(locAlpha0), curU0(:, locAlpha0), curRho, scaleFactor );
        %update for training data
        [ curZ0(:, locAlpha1) ] = updatez0_v4( 1, A.LINK_FUNC, curY(:, locAlpha1), curZ0(:, locAlpha1), D, ...
            curW(:, locAlpha1), curW0(locAlpha1), curU0(:, locAlpha1), curRho, scaleFactor, [], [] );
        %% update z1
        if verbose == 1
            fprintf( 'updating z1... ' );
        end
        preZ1 = curZ1;
        [ curZ1 ] = updatez1_v2( curZ1, curU1, curW, curRho, A.lambda );
        %% update z2
        if verbose == 1
            fprintf( 'updating z2... ' );
        end
        preZ2 = curZ2;
        [ curZ2 ] = updatez2_v1( curZ2, curU2, curRho, cRblk, curW, A.theta, A.L1Flag );
        %% update u0, u1
        if verbose == 1
            fprintf( 'updating u0, u1, u2...\n ' );
        end
        preY = D*curW(:,:) + repmat( curW0(:)', sLen, 1);
        %primal residual r0
        curR0 = curZ0 - preY;
        curU0 = curU0 + curRho*curR0;
        %primal residual r1
        curR1 = curZ1 - curW;
        curU1 = curU1 + curRho*curR1;
        %primal residual r2
        diffW = cRblk*curW';
        curR2 = curZ2 - diffW;
        curU2 = curU2 + curRho*curR2;
        %compute the primal residual all
        rf = norm( [curR0(:); curR1(:); curR2(:)] );
        pDim = length( curR0(:) ) + length( curR1(:) ) + length( curR2(:) );
        %compuate primal epsilon
        relEpsPri = max( norm( [curZ0(:); curZ1(:); curZ2(:) ]), norm( [preY(:); curW(:); diffW(:)] ) );
        epsPri = sqrt(pDim) * A.EPS_ABS_PRI + relEpsPri * A.EPS_REL_PRI;
        
        %compute the dual variable to see if it is converge
        %dual residual s0, s1, s2
        curS0 = curRho*( D'*( curZ0 - preZ0 ) );
        curS1 = curRho*( curZ1 - preZ1 );
        curS2 = curRho*( cRblk'*( curZ2 - preZ2 ) );
        sf = norm( curS0(:) + curS1(:) + curS2(:) );
        nDim = length( curS0(:) );
        %compuate dual epsilon
        tmp = D'*curU0;
        tmp = [tmp; sum(curU0)];
        if isempty( tmp )
            tmp = 0;
        end
        tmp2 = cRblk'*curU2;
        maxRel = max( [ norm( tmp(:) ), norm( u1(:) ), norm( tmp2(:) ) ] );
        epsDual = sqrt(nDim) * A.EPS_ABS_DUAL + maxRel * A.EPS_REL_DUAL;
        resRecAry(itNumADMM, :, j) = [rf, epsPri, sf, epsDual];
        changeW = full( max( abs( curW(:) -  preW(:) ) ) );
        if isempty(changeW)
            changeW = 0;
        end
        if isempty(curW)
            maxW = inf;
        else
            maxW = max(curW(:));
        end
        %         time = toc;
        %         fprintf( 'time: %g\n', time );
        %% update the next stage rho
%         if rf < epsPri && sf > epsDual
%             curRhoAry(itNumADMM+1, j) = 0.5* curRhoAry(itNumADMM, j);
%         elseif rf > epsPri && sf < epsDual
%              curRhoAry(itNumADMM+1, j) = 2* curRhoAry(itNumADMM, j);
        if mod(itNumADMM, 10 ) == 0
            curRhoAry(itNumADMM+1, j) = 10* curRhoAry(itNumADMM, j);
        elseif rf > 10*(epsPri/epsDual)*sf
            curRhoAry(itNumADMM+1, j) = 2* curRhoAry(itNumADMM, j);
        elseif sf > 10*(epsDual/epsPri)*rf
            curRhoAry(itNumADMM+1, j) = 0.5* curRhoAry(itNumADMM, j);    
        else
            curRhoAry(itNumADMM+1, j) = curRhoAry(itNumADMM, j);
        end
        
        %% return the variables
        A.W(usedMolIdx, loc) = curW;
        A.W0(loc) = curW0;
        A.z0(:, loc) = curZ0;
        u0(:, loc) = curU0;
        A.z1(usedMolIdx, loc) = curZ1;
        u1(:, loc) = curU1;
       
        A.z2{j}(:, usedMolIdx) = curZ2;
        u2{j} = curU2;
        [LPAryADMM(itNumADMM+1, j), dTerms, lTerms, ~, tTerms] = LP_DL_Poiss( A.LINK_FUNC, aMatrix, Y, A.W(usedMolIdx, :), A.W0, D, A.lambda, 0, A.theta, scaleFactor, logFY, 0 );
        
        %% breaking condition
        if  LPAryADMM(itNumADMM+1, j) < LPAryADMM(1, j) && ...
                ( rf < epsPri && sf < epsDual ) && ...
                ( changeW < maxW * A.W_TOL )
            %( max( abs(rf(:)) ) < 1e-3 && max( abs(sf(:)) ) < 1e-3 ) || ...
            break;
        end
        
        %correct the RHO_INIT with the ratio of data terms and penalty
        %terms in the first iteration.
        if itNumADMM == 1
            %set rho to be at least larger then 1e-4
            if (dTerms / (lTerms+tTerms) >= 0.5 && dTerms / (lTerms+tTerms) <= 10 ) && curRhoAry(1, j) < 1e-4
                curRhoAry(itNumADMM+1, j) = 1e-4;
            %set rho to be at least larger then 1e-6
            elseif (dTerms / (lTerms+tTerms) >= 50 && dTerms / (lTerms+tTerms) <= 1000 ) && curRhoAry(1, j) < 1e-6
                curRhoAry(itNumADMM+1, j) = 1e-6;
            end
        end
        if j == BlkDS.blkNum
            if A.WHistFlag == 1
                A.WHistCell{itNumADMM} = A.W;
            end
        end
        if verbose == 1
            fprintf( '(%d: %g) %g %g %g %g %g %g\n',itNumADMM, LPAryADMM(itNumADMM+1, j), rf, epsPri, sf, epsDual, curRhoAry(itNumADMM, j), changeW );
        end
        if itNumADMM == A.itNum
            if LPAryADMM(itNumADMM+1, j) >= LPAryADMM(1, j)
                if A.itNum < oITNum*2
                    A.itNum = A.itNum + 5;
                    tmp = LPAryADMM;
                    LPAryADMM = zeros( A.itNum+1, BlkDS.blkNum );
                    LPAryADMM(1:(end-5),:) = tmp;
                else
                    stuckFlag = 1;
                    break;
                end
            end
        end
        itNumADMM = itNumADMM + 1;
    end
    if stuckFlag == 1
        break;
    end
end
%return to original setting of itNum for next iteration DL.
A.itNum = oITNum;
A.LPAry = LPAryADMM;
A.rhoAry = curRhoAry;
A.resAry = resRecAry;
A.WDiff = changeW;
A.stuckFlag = stuckFlag;
end

function C = replicateC( inD, repNum )
[i, j] = find( inD ~= 0 );
[sLen, mLen] = size(inD);
valVec = nonzeros( inD );
valVec = repmat( valVec, repNum, 1 );
yVec = zeros(length(i)*repNum, 1);
xVec = zeros(length(i)*repNum, 1);
unitLen = length(i);
for it = 1:repNum
    xVec( ((it-1)*unitLen+1):(unitLen*it) ) = j + ((it-1)*mLen);
    yVec( ((it-1)*unitLen+1):(unitLen*it) ) = i + ((it-1)*sLen);
end
C = sparse( yVec, xVec, valVec, sLen*repNum, mLen*repNum );

% for it = 2:repNum
%     tmp = S( ((it-1)*sLen+1):(it*sLen), : ) - S(1:sLen, :);
%     if ~isempty( find( tmp ~= 0, 1 ) )
%         fprintf( '%d\n', it );
%     end
% end
end

function C = rowCombSparseMatrix( A, B )
[i, j] = find( A ~= 0 );
[sLenA, mLenA] = size(A);
[q, r] = find( B ~= 0 );
[sLenB, mLenB] = size(B);
if mLenA ~= mLenB
    fprintf( 'row combine sparse matrix error.\n' );
    C = [];
else
    valVecA = nonzeros( A );
    valVecB = nonzeros( B );
    allVec = [valVecA; valVecB];
    q = q + sLenA;
    allYAxis = [i; q];
    allXAxis = [j; r];
    C = sparse( allYAxis, allXAxis, allVec, sLenA+sLenB, mLenA );
end
end

function [C, target] = modifyCWithSparsity( CforUW, newWInfo, loc, mLen, nrPart1 )
for i = 1:length(loc)
    cLoc = loc(i);
    cW = newWInfo{cLoc};
    cW = cW+(i-1)*mLen;
    cWRange = (1+(i-1)*mLen):(i*mLen);
    target = setdiff(cWRange, cW);
    CforUW(1:nrPart1, target) = 0;
end
[a, b] = find( CforUW((nrPart1+1):size(CforUW, 1), :)' ~= 0 );
target = [];
for i = 1:2:length(a)
    [q, r] = ind2sub([mLen, length(loc)], a(i:(i+1)) );
    if isempty( find( newWInfo{loc(r(1))} == q(1), 1 ) ) || isempty( find( newWInfo{loc(r(2))} == q(2), 1 ) )
        target = [ target b(i)];
    end
end
target = target + nrPart1;
CforUW( target, :) = 0;
C = CforUW;
end
