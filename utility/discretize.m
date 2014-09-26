function [ da_vec, centers ] = discretize( a_vec, binNum )
%discretize discretize a vector according to user's input bin number
[nelements,centers] = hist( a_vec, binNum );
da_vec = zeros( size(a_vec) );
for i = 1:length(a_vec)
    [~,idx] = min( abs(a_vec(i)-centers) );
    da_vec(i) = centers(idx);
end

