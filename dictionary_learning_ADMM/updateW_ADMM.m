function [WResStruct] = updateW_ADMM( Y, D, W, W0, z0, z1, aMatrix, BlkDS, itNum, lambda, phi, theta, scaleFactor, logFY )

%% variable setting
fprintf( 'iteration number for W_update_ADMM = %d\n', itNum );
temp = (aMatrix==0)&(BlkDS.indMap==1);
fprintf( 'held-out sample number = %d, training sample number = %d\n', length( find( temp == 1 ) ), length( find( aMatrix == 1 ) ) );
fprintf( 'lambda = %d, phi = %d, theta = %d\n', lambda, phi, theta );
z0 = z0(:, :); z1 = z1(:, :); W = W(:, :); W0 = W0(:, :);
LPAryADMM = zeros( itNum+1, 1 );
resRecAry = zeros( itNum, 4 );
curRhoAry = zeros( itNum, 1 );
curRhoAry(1) = 1;
u0 = sparse( zeros( size( z0 ) ) );
u1 = sparse( zeros( size( z1 ) ) );
LPAryADMM(1) = LP_DL_Poiss( aMatrix, Y, W, W0, D, lambda, phi, theta, scaleFactor, logFY, 0 );
[sLen, hei, wid] = size( Y );
nLen = hei*wid;
[~, mLen] = size( D );
indMap = BlkDS.indMap;

%% contructing genSparseGroupingMatrix
[Rall] = genSparseGroupingMatrix( hei, wid, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end

%% main optimization process
for itNumADMM = 1:itNum
    prevW_ADMM = W;
    fprintf( 'LP: %g ', LPAryADMM(itNumADMM) );
    curRho = curRhoAry(itNumADMM);
%% update W
    fprintf( 'updating W...\t' );
    moutD = [D ones(sLen, 1)];
    tmpEye = eye(mLen+1, mLen+1);
    tmpEye(end, end) = 0;
    CforUW = sparse( [moutD; tmpEye] );
    tic
    parfor i = 1:nLen
        if  indMap(i)
            [ w_w0 ] = update_w_w0_v2( z0(:,i), u0(:,i), z1(:,i), u1(:,i), CforUW, curRho, [W(:,i); W0(i)] );
            W(:,i) = w_w0(1:end-1); W0(i) = w_w0(end);
        end
    end

%% testing the algorithm combining all samples as a matrix to optimize the weight. currently this function is disable.
%     %construct Big CforUW
%     [i, j] = find( CforUW ~= 0 );
%     validSam = length( find( indMap==1 ) );
%     iVec = zeros( length(i)*validSam, 1 );
%     jVec = zeros( length(j)*validSam, 1 );
%     valVec = zeros( length(j)*validSam, 1 );
%     for z = 1:validSam
%         iVec( ((z-1)*length(i)+1):z*length(i) ) = i;
%         jVec( ((z-1)*length(j)+1):z*length(j) ) = j;
%         valVec( ((z-1)*length(j)+1):z*length(j) ) =  CforUW(CforUW>0);
%     end
%     BCforUW = sparse( iVec, jVec, valVec, size(CforUW, 1)*validSam, size(CforUW, 2)*validSam );
%     %construct Big W_init, Big Dictionary
%     BW_init = zeros( (mLen+1)*validSam, 1 );
%     BD = zeros( (size(z0, 1)+size(z1, 1)+1)*validSam, 1 );
%     z = 1;
%     for i = 1:nLen
%         if indMap(i)
%             BW_init( ((z-1)*(mLen+1)+1):z*(mLen+1) ) = [W(:,i);W0(i)];
%             BD( ((z-1)* (size(z0, 1)+size(z1, 1)+1)+1):z*(size(z0, 1)+size(z1, 1)+1)) = [z0(:, i)+1/curRho*u0(:, i); z1(:, i)+1/curRho*u1(:, i); 0];
%             z = z + 1;
%         end
%     end
%     BD = sparse(BD);
%     
%     %use lsqlin to optimze, the process is slow
%     options = optimoptions('lsqlin');
%     options.MaxIter = 100;
%     options.TolFun=1e-6;
%     options.Display='none';
%     lb = zeros(size(BCforUW, 2), 1);
%     %     tic
%     % [w_w0,~,~,~, ~] = lsqlin( BCforUW, BD, [], [], [], [], lb, [], BW_init, options );
%     % toc;
%     
%     %yse nnLeastR to optimize, the process is fast
%     opts=[];
%     opts.mFlag=0;       % treating it as compositive function
%     opts.lFlag=0;       % Nemirovski's line search
%     opts.init=2;        % starting from a zero point
%     % termination criterion
%     opts.tFlag=5;       % run .maxIter iterations
%     opts.maxIter=100;   % maximum number of iterations
%     % normalization
%     opts.nFlag=0;       % without normalization
%     % regularization
%     % opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
%     tic;
%     [x2, funVal2]= nnLeastR( BCforUW, BD,1e-6, opts);
%     toc;
%% update z0
    fprintf( 'updating z0...\t' );
    prevZ0 = z0;    %for computing the dual residual
    parfor i = 1:nLen
        if indMap(i)
            [ z0(:, i) ] = updatez0_v3( aMatrix(i), Y(:,i), z0(:,i), D, W(:, i), W0(i), u0(:,i), curRho, scaleFactor );
        end
    end

%% update z1
    fprintf( 'updating z1...\n' );
    prevZ1 = z1;    %for computing the dual residual
    for j = 1:BlkDS.blkNum
        loc = BlkDS.B2GMap{j};
        curW = W(:, loc);
        curU1 = u1(:, loc);
        tmpZ1 = zeros( size(W, 1), length( loc ) );
        curR = Rblk{j};
        parfor i = 1:size(W, 1)
            %[ tmpZm ] = updatez1_m( [], [], curRho, curW(i, :), curU1(i, :), lambda+1e-8, theta+1e-8, curR, []);
            tmpZm = updatez1_m_ridgePen_SLEP( curRho, curW(i, :), curU1(i, :), lambda, theta, curR );
            tmpZ1(i, :) = tmpZm;
        end
        z1(:, loc) = tmpZ1;
    end

%% update u0, u1
    preY = D*W(:,:) + repmat( W0(:)', sLen, 1);
    %primal residual r1
    r1 = z0 - preY;
    %         [val, o] = max(r1(:));
    %         fprintf( '%g %g %g\n', val,  o, full(z0(o)));
    %         for i = 1:wid*hei
    %             if ~indMap(i) && ~isempty( find( r1(:,i) ~= 0, 1 ) )
    %                 fprintf( 'loc: %d', i );
    %                 r1(:, i) = 0;
    %             end
    %         end
    u0 = u0 + curRho*r1;
    %compute the primal residual all
    rf = r1;
    r2 = z1 - W;
    %         for i = 1:wid*hei
    %             if ~indMap(i) && ~isempty( find( r2(:,i) ~= 0, 1 ) )
    %                fprintf( 'loc: %d', i );
    %                 r2(:, i) = 0;
    %             end
    %         end
    %
    %         fprintf( 'r2 < 0 : %d\n', length(find(r2b(:)<0)) );
    u1 = u1 + curRho* r2;
    rf = [rf; r2];
    rfn2 = norm(rf(:), 2);

    %compute the dual variable to see if it is converge
    % %         dualVar = [u0; u1{1}];
    % %         dualVar = norm(dualVar(:), 2);
    %
    %dual residual s0
    s0 = curRho*( D'*( z0 - prevZ0 ) );
    % %         for i = 1:wid*hei
    % %             if ~indMap(i) && ~isempty( find( s0(:,i) ~= 0, 1 ) )
    % %                 fprintf( 'loc: %d', i );
    % %                 s0(:,i) = 0;
    % %             end
    % %         end
    s0 = [s0; sum( z0 -prevZ0, 1 )];
    s0n2 = norm( s0(:), 2 );
%
%         %compute the more dual residual values needed at the block condtion
% %         duResMore = [];
% %         for i = 1:hei*wid
% %             bID = BlkDS.mapW2B{i};
% %             if length(bID) > 1
% %                 tmp = 0;
% %                 for j = 2:length(bID)
% %                     tmp = tmp + curRho*( ( z1{bID(j)}(:, i)- prevZ1{bID(j)}(:, i) ) );
% %                 end
% %                 duResMore = [duResMore; tmp];
% %             end
% %         end
% %         duResMore = norm(duResMore, 2);
% 
%         %converge condition
%         %primal residual condition
    EPS_ABS = 1e-4;
    EPS_REL = 1e-3;
%         
    epsPri = sqrt( size(rf, 1)*size(rf, 2) )*EPS_ABS;
    relEpsPri = max( [norm( z0(:), 2 ) norm( preY(:) ) norm( z1(:), 2 ), norm( W(:), 2 )] );
    epsPri = epsPri + relEpsPri*EPS_REL;
%         
        %dual residual condition
%         for i = 1:BlkDS.blkNum
%             if dataMask(i) == 0
%                 continue;
%             else
%                 tmp = tmp + u1{i};
%             end
%         end
    tmp = D'*u0;
    tmp = [tmp; sum(u0)];
    maxRel = max( norm( tmp(:), 2 ), norm( u1(:), 2 ) );
    epsDual = sqrt( size(s0, 1) * size(s0, 2) )*EPS_ABS + maxRel*EPS_REL;
    
    resRecAry(itNumADMM, :) = [rfn2, epsPri, s0n2, epsDual];
    tmp = max( abs( W(:) -  prevW_ADMM(:) ) );
    fprintf( '%d: %g %g %g %g %g %g\n',itNumADMM, rfn2, epsPri, s0n2, epsDual, curRhoAry(itNumADMM), full(tmp) );
    LPAryADMM(itNumADMM+1) = LP_DL_Poiss( aMatrix, Y, W, W0, D, lambda, phi, theta, scaleFactor, logFY, 0 );
    
    if ( rfn2 < epsPri && s0n2 < epsDual ) ||...
            ( max( abs(rf(:)) ) < 1e-3 && max( abs(s0(:)) ) < 1e-3 ) || ...
            ( tmp < 1e-3 )
        break;
    end

    %update the next stage rho
    if rfn2 > 10*s0n2
        curRhoAry(itNumADMM+1) = 2 * curRhoAry(itNumADMM);
        %         elseif s0n2 > 2*rfn2
        %             curRhoAry(itNumADMM+1) = 0.677 * curRhoAry(itNumADMM);
    elseif s0n2 > 10*rfn2
        curRhoAry(itNumADMM+1) = 0.5 * curRhoAry(itNumADMM);
    else
        curRhoAry(itNumADMM+1) = curRhoAry(itNumADMM);
    end
    toc;
%     rhoCell{it} = curRhoAry;
%     resRecCell{it} = resRecAry;
end
WResStruct.W = sparse( W );
WResStruct.W0 = sparse( W0 );
WResStruct.z0 = sparse( z0 );
WResStruct.z1 = sparse( z1 );
WResStruct.LPAry = sparse( LPAryADMM );
WResStruct.rhoAry = sparse( curRhoAry );
WResStruct.resAry = sparse( resRecAry );
WResStruct.WDiff = sparse( full(tmp) );
end
