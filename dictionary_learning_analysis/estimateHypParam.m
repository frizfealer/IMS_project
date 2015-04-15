function [ dfRang, pRang ] = estimateHypParam( varFilePath, mode, varargin )
%--------------------------------------------------------------------------
% estimateHypParam: estimate appropriate hyper parameter ranges
%--------------------------------------------------------------------------
% DESCRIPTION:
%   estimate appropriate hyper parameter ranges by checking 
%	the resulting degree of freedom
%
% INPUT ARGUMENTS:
%   varFilePath, basic variable path
%	mode, 1 for lambda, 2 for phi, and 3 for theta
%	varargin{1} varargin{2}, 
%		the smallest and the largest value of your hyperparameters
%		(in log10 domain )
% OUTPUT ARGUMENTS:
%   pRang, hyperparameter range
%	dfRang, degree of freedom range
load( varFilePath );
[ param ] = setDLParameters();
param.OUTER_IT_NUM = 1;
%% testing lambda
if mode == 1
	if ~isempty( varargin )
		lmin = varargin{1}; lmax = varargin{2};
	else
		lmin = -1; lmax = 4;
	end
	lVec = logspace( lmin, lmax, 10);
    lVec = [0 lVec];
    dfVec1 = ones( length(lVec), 1 )*-1;
    for i = 1:length(lVec)
        param.UP_D_IT_NUM = 1;
        [ expRec ] = dictionaryLearning_ADMM_v6_2( IMSD.dataCube, [], DTemplate2, [], lVec(i), 1e-32, 1e-32, aMatrix, IMSD.BlkDS, [], 'temp.mat', param );
        if expRec.WStuckFlag == 0
            dfVec1(i) = length(find(expRec.W(:)>1e-2));
        end
    end
    figure; plot(dfVec1); title( 'lambda v.s. d.f.' );
    dfRang = dfVec1;
	pRang = lVec;
%% testing phi
elseif mode == 2
	if ~isempty( varargin )
		pmin = varargin{1}; pmax = varargin{2};
	else
		pmin = 1; pmax = 10;
	end
    pVec = logspace( pmin, pmax, 10);
    pVec = [0 pVec];
    dfVec2 = ones( length(pVec), 1 )*-1;
    for i = 1:length(pVec)
        param.UP_D_IT_NUM = 200;
        [ expRec ] = dictionaryLearning_ADMM_v6_2( IMSD.dataCube, [], DTemplate2, [], 1e-32, 1e-32, pVec(i), aMatrix, IMSD.BlkDS, [], 'temp.mat', param );
        if expRec.DStuckFlag == 0
            dfVec2(i) = length(find(expRec.D(:)>1e-2));
        end
    end
    figure; plot(dfVec2); title( 'phi v.s. d.f.' );
    dfRang = dfVec2;
	pRang = pVec;
    %% testing theta
elseif mode == 3
    [~, hei, wid] = size(IMSD.dataCube);
    [Rall] = genSparseGroupingMatrix( hei, wid, 1 );
    for i = 1:IMSD.BlkDS.blkNum
        Rblk{i} = [];
        Rb = Rall(:, IMSD.BlkDS.B2GMap{i});
        ins = sum(abs(Rb), 2);
        Rb(ins<2,:) = [];
        Rblk{i} = Rb;
    end
	if ~isempty( varargin )
		tmin = varargin{1}; tmax = varargin{2};
	else
		tmin = -4; tmax = -2;
	end
    tVec = logspace( tmin, tmax, 10 );
    tVec = [0 tVec];
    dfVec3 = ones( length(tVec), 1 )*-1;
    for i = 1:length(tVec)
        param.UP_D_IT_NUM = 1;
        [ expRec ] = dictionaryLearning_ADMM_v6_2( IMSD.dataCube, [], DTemplate2, [], 1e-32, tVec(i), 1e-32, aMatrix, IMSD.BlkDS, [], 'temp.mat', param );
        if expRec.WStuckFlag == 0
            for j = 1:IMSD.BlkDS.blkNum
                loc = IMSD.BlkDS.B2GMap{j};
                tmp = expRec.W(:, loc);
                tmp = Rblk{j}*tmp';
                dfVec3(i) = dfVec3(i) + norm( tmp(:) )^2;
            end
        end
        figure; imagesc2( reshape( max( expRec.W, [], 1 ), hei, wid ) );
    end
    dfRang = dfVec3;
	pRang = tVec;
end
end

