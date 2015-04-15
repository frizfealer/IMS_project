function [index] = findMatchedPeaks( mzAxis, targetMZ, tol )
%--------------------------------------------------------------------------
% findMathcedPeaks: find the index of a m/z values in the data
%--------------------------------------------------------------------------
% DESCRIPTION:
%   a short function to return the index in the m/z list of the data with a
%   target m/z value.
%
% INPUT ARGUMENTS:
%   mzAxis, a vector of m/z values in the data
%   targetMZ, a target m/z value we want to find.
%   tol, tolerance, the default setting is 500 ppm.
% OUTPUT ARGUMENTS:
%   index, if found, it returns an index in mzAxis; otherwise returns -1.
[val, idx] = min( abs(targetMZ - mzAxis ) );
if val <= targetMZ*tol
    index = idx;
else
    index = -1;
end
end