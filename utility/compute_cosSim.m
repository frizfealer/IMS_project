function [ val ] = compute_cosSim( a, b )
%compute_cosSim compute cosine similarity
normA = norm( a );
normB = norm( b );
val = sum(a.*b)/(normA*normB);
end

