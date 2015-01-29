function [ gY ] = genY_differentNoise( LINK_FUNC, gD, gW, gW0, varargin )
%genY generate Y from gD, gW, gW0 with poisson noise
%   Detailed explanation goes here
[SLEN] = size(gD, 1);
[~,HEIGHT, WIDTH] = size( gW );
gY = gD*gW(:,:) + repmat( gW0(:)', SLEN, 1 );
gY = reshape( gY, SLEN, HEIGHT, WIDTH );
% gY = gY + rand( sLen, IHEIGHT, IWIDTH );
if strcmp( LINK_FUNC, 'log' ) == 1
    lambda = exp( gY );
    gY = poissrnd( lambda );
    %if no weight on that location, let the signal be zeros
    ins = sum( gW, 1 );
    gY(:,ins==0) = 0;
elseif strcmp( LINK_FUNC, 'identity' ) == 1
    lambda = gY;
    gY = poissrnd( lambda );
    %if no weight on that location, let the signal be zeros
    ins = sum( gW, 1 );
    gY(:,ins==0) = 0;
elseif strcmp( LINK_FUNC, 'log_gaussain' ) == 1
	lambda = gY;
    gY = exp( randn( size(lambda) )/1e3 + lambda );
elseif strcmp( LINK_FUNC, 'negative_binomial' ) == 1
    lambda = exp( gY );
    %R: 1/alpha, P: 1/(1+alpha*lambda)
    if length(varargin) == 1
        kappa = varargin{1};
    else
        kappa = 1e-2;
    end
    R = 1/kappa;
    P = 1/(1+kappa*lambda);
    gY = nbinrnd(R,P);
    ins = sum( gW, 1 );
    gY(:,ins==0) = 0;
end
if ~isempty( find( isinf(lambda), 1 ) )
    fprintf( 'the value of w is too high.\n' );
end

end

