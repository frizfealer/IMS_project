function [ resVec ] = evaluateRecovery( dict, goodPeaks, badPeaks, thres )
%evaluateRecovery evaluate each compounds in the dict.
%return a score for each compounds, that if it has a good peaks above a
%threshold, add 1, if it has a bad peaks above a threshold, decrease 1.
%goodPeaks, badPeaks, the number of peaks in dict

resVec = zeros( size(dict, 2), 1);
assert( isempty( setdiff( goodPeaks, 1:size(dict, 1) ) ) );
assert( isempty( setdiff( badPeaks, 1:size(dict, 1) ) ) );
for i = 1:size(dict, 2)
    idx = find( dict(:, i) >= thres );
    resVec(i) = resVec(i) + length( intersect( goodPeaks, idx ) );
    resVec(i) = resVec(i) - length( intersect( badPeaks, idx ) );
end


end

