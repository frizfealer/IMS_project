function [ mVec ] = computeMapping( vec1, vec2, eps )
%computeMapping Compute every value in vec1 to the nearest values in vec2
%vec1, vec2: vectors
%eps the error range expected when mapping
%return mVec, a two dimension vector [n,2], such that
%[n,1] is the value in vec2, and [n, 2] is the index in vec2
mVec = repmat( vec1, 1, 2 );
for i = 1:length(vec1)
    [ mVal, idx] = min( abs( vec1(i) - vec2 ) );
    if mVal > eps
        mVec(i, :) = -nan;
    else
        mVec(i, 1) = vec2(idx);
        mVec(i, 2) = idx;
    end
end

end

