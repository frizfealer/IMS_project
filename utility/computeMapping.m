function [ mVec ] = computeMapping( vec1, vec2, eps )
%--------------------------------------------------------------------------
% computeMapping: compute a mapping from vec1 to vec2, with eps ppm.
%--------------------------------------------------------------------------
%
% INPUT ARGUMENTS:
%   vec1, the target vector you want to compute mapping, usually has a size
%   smaller than vec2.
%   vec2, the reference vector you want to compute mapping.
%   eps, the error tolerance in ppm sense that related to the target vector
%   value.
% OUTPUT ARGUMENTS:
%   mVec, a two dimension vector [n,2], such that [n,1] is the value in 
%   vec2, and [n, 2] is the index in vec2 (n is the length of vec1). If
%   there is no mapping for a point in vec1, the corresponding row in mVec
%   is -nan.
mVec = repmat( vec1, 1, 2 );
for i = 1:length(vec1)
    [ mVal, idx] = min( abs( vec1(i) - vec2 ) );
    if mVal > vec1(i)*eps
        mVec(i, :) = -nan;
    else
        mVec(i, 1) = vec2(idx);
        mVec(i, 2) = idx;
    end
end

end

