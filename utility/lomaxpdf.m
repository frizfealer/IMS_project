function [ val ] = lomaxpdf( alpha, lambda, x )
%lomaxpdf pdf of lomax probability distribution
% alpha/lambda*(1+x/lambda)^-(alpha+1)
val = alpha/lambda * (1+x./lambda).^-(alpha+1);

end

