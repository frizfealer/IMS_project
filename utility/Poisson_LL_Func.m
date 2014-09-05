function [ val ] = Poisson_LL_Func( Y, D, W, preY )
%Poisson_LL_Func compute Poisson log likelihood function
%indMap, a index map for samples, if respected entry is 0, it means the
%sample is not used in log likelihood computation
if isempty( preY )
    preY = D*W;
end
val = Y.*preY - exp( preY );
val = sum( sum( val ) );
end

