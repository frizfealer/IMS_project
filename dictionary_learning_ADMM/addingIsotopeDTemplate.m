function [ uDTemplate, uDIonName ] = addingIsotopeDTemplate( mzAxis, errRange, DTemplate, DIonName, ISONumber, usedIonForms )
%--------------------------------------------------------------------------
% addingIsotopeDTemplate: adding Isotope ions to the dictionary template
%--------------------------------------------------------------------------
% DESCRIPTION:
%   Adding more possible ion types that comes from isotopes of M to the
%   dictionary template.
%
% INPUT ARGUMENTS:
%   mzAxis, the mzAxis of the dictionary.
%   errRange, the error tolerance of m/z shifting, should be consistenet
%   with binning error. Default setting is 5e-4 (500 ppm).
%   DTemplate, DIonName, SpeciesM, the output a dictioanry template 
%   generated from genDTemplate
%   ISONumber, is the number of isotope to use, shoud be 1 or 2 ([M+1],
%   [M+2])
%   usedIonForms, specify what ion forms are allowed to add isotype.
%   E.g. [1 3 5] means adding isotypes only on the first, the third, and
%   the fifth ion types. CURRENTLY DEPRECATED.
% OUTPUT ARGUMENTS:
%   DTemplate, a template with 1 means possible and 0 means impossible
%   entries of a dictionary.
%   DIonName, all possible adducts generated that signal. Every dictionary
%   elements has a cell of possilbe DIonNames and each entry in this cell
%   correspondes to one possible set of adducts.
%   speciesM, all possible molecular weight generate that signal. Every
%   dictionary elements has a vector of possible molecular weights.

if ISONumber == 1
    mzDiff = [1];
elseif ISONumber == 2
    mzDiff = [1 2]';
end

uDTemplate = DTemplate;
uDIonName = DIonName;
for i = 1:size(DTemplate, 2)
    idx = find( DTemplate(:, i) == 1 );
    cName = DIonName{i};
    mapTab = cell( 1, 2 );
    cnt = 1;
    for j = 1:length(idx)
        mapTab(cnt, :) = [idx(j), cName(j)]; cnt = cnt + 1;
        cMZ = mzAxis(idx(j));
        possibleRelatedPeaks = cMZ + mzDiff;
        [ mVec ] = computeMapping( possibleRelatedPeaks, mzAxis, errRange );
        curPeak = mVec( :, 2 );
        if length(mzDiff) == 2
            if isnan(curPeak(1)) && ~isnan(curPeak(2))
                curPeak(2) = nan;
            end
        end
        for z = 1:length(curPeak)
            if ~isnan(curPeak(z))
                mapTab(cnt, :) = [curPeak(z), {[cName{j} '+' num2str(z)]}];
                cnt = cnt + 1;
            end
        end
    end
    [~, idx] = sort( cell2mat(mapTab(:,1)) );
    uDTemplate(cell2mat(mapTab(:,1)), i) = 1;
    tmp = mapTab(:, 2);
    uDIonName{i} = tmp(idx);
end

end