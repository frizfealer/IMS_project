function C = fastIntersect( A, B )
%fastIntersect Intersect of two sets A, B
%   A and B must be non-negative Integers 
% and A and B must be unique.
if ~isempty(A)&&~isempty(B)
   P = zeros(1, max(max(A),max(B)) ) ;
   P(A) = 1;
   C = B(logical(P(B)));
else
    C = [];

end