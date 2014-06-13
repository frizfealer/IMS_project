function [ rD ] = updateD_v4_cvx( inY, outW, outW0, D_init, sLen, mLen, DTemplate, aMatrix, lambda, theta, scaleFactor, logFactorialY )
%updateD_v4_cvx update the whole dictionary using fmincon without ADMM
%logFactorialY is the log factorial of input Y, the same size of Y
%scaleFactor set to 1.
assert( sLen == size(D_init, 1) );
assert( mLen == size(D_init, 2) );


Y = inY(:, :); %Y size[sLen h*w] 
W = outW(:, :); %W size[mLen h*w]
prevD = D_init; %D size[sLen mLen]
preLP = LP_DL_Poiss( aMatrix, inY, outW, outW0, D_init, lambda, 0, theta, scaleFactor, logFactorialY );
threshold = 1e-10*length(find(DTemplate(:)~=0));

staticTerm = repmat( outW0(:)', sLen, 1 );
zeroEnt = ( DTemplate == 0 ); %get the indices of entries in dictioanry that should be zero
nonZeroEnt = (DTemplate == 1 );

%% cvx code start
cvx_begin
% cvx_solver sedumi
variable uD(sLen, mLen)
% expression len
% len = sum(uD.^2,1).^0.5;
minimize( sum( sum ( -Y.*(staticTerm+uD*W) + exp((staticTerm+uD*W)) +  logFactorialY) ) )
subject to
uD(zeroEnt) == 0;
uD(nonZeroEnt) >= 0;
for i = 1:mLen
    norm( uD(:, i ) ) <= 1;
end
cvx_end


rD = uD;
    curLP = LP_DL_Poiss( aMatrix, inY, outW, outW0, rD, lambda, 0, theta, scaleFactor, logFactorialY );
    fprintf( '%g %g %g\n', norm( rD(:) - prevD(:), 2 ), threshold, curLP-preLP );
end

