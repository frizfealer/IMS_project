function val = logfactorial_e( inY, MAXNUM )
%inY a vectors
%return val the same size as inY
val = zeros( size( inY ) );
maxVal = max( inY(:) );
if maxVal > MAXNUM %about using  1e8, 0.7451 GB
    maxVal = 1e8;
end
maxAry = zeros( maxVal, 1);
for i = 1:maxVal
    maxAry(i) = log(i);
end;
for i = 1:size(inY(:))
    if inY(i) <= maxVal
        val(i) = sum( maxAry(1:inY(i)) );
    else
        val(i) = sum( maxAry(1:end) );
        for j = (maxVal+1):inY(i)
            val(i) = val(i) + log(j);
        end
    end
end
end
