function [ tIdxVec ] = check_related_peaks_index( mzTarget, mzAxis, W )
%check_related_peaks_index check related peaks index in mzAxis, with
%TOLERANCE 1Da
TOLERANCE = 1;
tIdxVec = zeros( length(mzTarget), 1 );
for i = 1:length(mzTarget)
    [~,t1] = min(abs(mzAxis-mzTarget(i)));
    insVec = [];
    idxVec = [];
    insVec(1) = max(W(t1, :));
    idxVec(1) = t1;
    if t1 - 1 > 0 && abs(mzAxis( t1 - 1 ) - mzTarget(i)) < TOLERANCE
        insVec(end+1) = max(W(t1-1,:));
        idxVec(end+1) = t1-1;
    end
    if t1 + 1 <= length(mzAxis) && abs( mzAxis( t1+1 ) - mzTarget(i) ) < TOLERANCE
        insVec(end+1) = max(W(t1+1,:));
        idxVec(end+1) = t1+1;
    end
    [~,ii] = max(insVec);
    tIdxVec(i) = idxVec(ii);
end
end

