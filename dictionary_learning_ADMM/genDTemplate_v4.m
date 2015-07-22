function [ DTemplate, DIonName, speciesM ] = genDTemplate_v4( MPPATH, IMSD, errRange )
%--------------------------------------------------------------------------
% genDTemplate_v4: generate DTemplate version 4
%--------------------------------------------------------------------------
% DESCRIPTION:
%   A new way to generate dictionary templates. It allows only the local
%   maximun m/z has all possible ion types. Other m/z values can only have
%   M+H ion types. This function could be changed soon.
%
% INPUT ARGUMENTS:
%   MPPATH: the file path of ion types
%   errRange, the error tolerance of m/z shifting, should be consistenet
%   with binning error. Default setting is 5e-4 (500 ppm).
%   IMSD: the input data cube types.
% OUTPUT ARGUMENTS:
%   DTemplate, a template with 1 means possible and 0 means impossible
%   entries of a dictionary.
%   DIonName, all possible adducts generated that signal. Every dictionary
%   elements has a cell of possilbe DIonNames and each entry in this cell
%   correspondes to one possible set of adducts.
%   speciesM, all possible molecular weight generate that signal. Every
%   dictionary elements has a vector of possible molecular weights.

ins = max( IMSD.dataCube(:, :), [], 2 );
[~,locs] = findpeaks(ins);
[ DTemplate3, DIonName3, speciesM3 ] = genDTemplate_v3( IMSD.mzAxis, MPPATH, errRange, [], [], IMSD.mzAxis(locs)  );
[ DTemplate3, DIonName3, SpeciesM3 ] = improveDTemplate( 'positive', DTemplate3, DIonName3, speciesM3 );
ins = sum( DTemplate3, 2);
idx =  ins == 0 ;
[ DTemplate4, DIonName4, speciesM4 ] = genDTemplate_v3( IMSD.mzAxis, MPPATH, errRange, [1], [], IMSD.mzAxis(idx)  );
[ DTemplate4, DIonName4, SpeciesM4 ] = improveDTemplate( 'positive', DTemplate4, DIonName4, speciesM4 );
DTemplate = [DTemplate3 DTemplate4];
DIonName = [DIonName3, DIonName4];
speciesM = [SpeciesM3, SpeciesM4]; 

mvVec = zeros( length(speciesM), 1 );
for i = 1:length(speciesM)
    mvVec(i) = mean(speciesM{i});
end
[~, idx] = sort( mvVec );
speciesM = speciesM(idx);
DIonName = DIonName(idx);
DTemplate = DTemplate(:, idx);


end