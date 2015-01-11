function [ gY ] = genY_Poisson( LINK_FUNC, gD, gW, gW0 )
%genY generate Y from gD, gW, gW0 with poisson noise
%   Detailed explanation goes here
[SLEN] = size(gD, 1);
[~,HEIGHT, WIDTH] = size( gW );
gY = gD*gW(:,:) + repmat( gW0(:)', SLEN, 1 );
gY = reshape( gY, SLEN, HEIGHT, WIDTH );
% gY = gY + rand( sLen, IHEIGHT, IWIDTH );
if strcmp( LINK_FUNC, 'log' ) == 1
    lambda = exp( gY );
elseif strcmp( LINK_FUNC, 'identity' ) == 1
    lambda = gY;
end
if ~isempty( find( isinf(lambda), 1 ) )
    fprintf( 'the value of w is too high.\n' );
end
gY = poissrnd( lambda );
%if no weight on that location, let the signal be zeros
ins = sum( gW, 1 );
gY(:,ins==0) = 0;

end

