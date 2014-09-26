function [ val ] = computeMI( vec1, vec2, binVal )
%computeMI compute mutual information of vec1 and vec2
%with binVal is a vector of binning center values of both vectors
[nelement1,~] = hist( vec1, binVal );
[nelement2,~] = hist( vec2, binVal );
s1 = sum(nelement1);
s2 = sum(nelement2);
p12 = zeros( 1, 3 );
p1 = nelement1 / s1;
p2 = nelement2 / s2;

cLen = 1;
for i = 1:length(vec1)
    idx1 = ( vec1(i) == p12(:, 1) );
    idx2 = ( vec2(i) == p12(:, 2) );
    idx = idx1 & idx2;
    if ~isempty(find(idx==1, 1))
        p12(idx==1,3) = p12(idx==1, 3) + 1;
    else
        p12(cLen, 1) = vec1(i);
        p12(cLen, 2) = vec2(i);
        p12(cLen, 3) = p12(cLen, 3) + 1;
        cLen = cLen + 1;
    end
end
p12(:,3) = p12(:,3) / length(binVal)^2;

val = 0;
for i = 1:size( p12, 1 )
    bIdx1 =  p12(i, 1) == binVal ;
    bIdx2 =  p12(i, 2) == binVal ;
    val = val + p12(i, 3)*( log( p12(i,3) ) - log( p1(bIdx1) ) - log(p2(bIdx2) ) );
end
   
end

