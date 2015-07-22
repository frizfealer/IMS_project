function [ DTemplate, DIonName, speciesM ] = genDTemplate_v3( inpeakMZ, MPPATH, errRange, varargin )
%--------------------------------------------------------------------------
% genDTemplate_v3: generate DTemplate version 3
%--------------------------------------------------------------------------
% DESCRIPTION:
%   Generate dictionary Template with given ion types. A dictionary
%   template is a matrix of the same size of D with its sparsity pattern
%   encodes.
%
% INPUT ARGUMENTS:
%   inpeakMZ, [sLen 1], the m/z values used to generate D Template, should
%   be the input m/z values
%   MPPATH: the file path of ion types
%   errRange, the error tolerance of m/z shifting, should be consistenet
%   with binning error. Default setting is 5e-4 (500 ppm).
%   usedIonForms, an additiona parameter that users specify to generate m/z
%   differences. E.g. [1 3 5] means using only the first, the third, and
%   the fifth ion types.
%   removeRedFlag, a flag indicates wheter to remove redundant dictionary
%   elements or not. A redundant between two dictionary elements is as
%   follows. 
%   A:[1 1 0 0]
%   B:[1 1 1 0]
%   A is redudent of B.
%   The default setting is off. This flag can be on if speed is too slow.
% OUTPUT ARGUMENTS:
%   DTemplate, a template with 1 means possible and 0 means impossible
%   entries of a dictionary.
%   DIonName, all possible adducts generated that signal. Every dictionary
%   elements has a cell of possilbe DIonNames and each entry in this cell
%   correspondes to one possible set of adducts.
%   speciesM, all possible molecular weight generate that signal. Every
%   dictionary elements has a vector of possible molecular weights.

fileID = fopen(MPPATH);
C = textscan(fileID, '%f %s', 'delimiter', {','});
fclose(fileID);
mpAry = C{1,1};
nameAry = C{1,2};
[mpAry, idx] = sort( mpAry );
nameAry = nameAry(idx);

SLEN = length(inpeakMZ);
if length(varargin) >= 1 && ~isempty(varargin{1})
    usedIonForms = varargin{1};
else
    usedIonForms = 1:length(mpAry);
end
if length(varargin) >= 2 && ~isempty(varargin{2})
    removeRedFlag = varargin{2};
else
    removeRedFlag = 0;
end
if length(varargin) == 3 && ~isempty(varargin{3})
    targetPeakMZ = varargin{3};
else
    targetPeakMZ = inpeakMZ;
end
if isempty( errRange )
    errRange = 5e-4;
end

assert( max(usedIonForms) <= length( mpAry ) );
assert( min(usedIonForms) >= 1 );
interval = mpAry - mpAry(1);
mapping = zeros( length(usedIonForms), length(interval) );
for i = 1:length(usedIonForms)
    mapping(i,:) = interval - interval(usedIonForms(i));
end

speciesM = [];
DTemplate = [];
DIonName = {};
cnt = 1;
fprintf( 'generating DTemplate...\n' );
for i = 1:length(targetPeakMZ)
    %     fprintf( 'm/z %g\n', mMVec(i));
    curMZ = targetPeakMZ(i);
    for j = 1:length(usedIonForms)
        cMapping = mapping(j, :);
        possibleRelatedPeaks = cMapping+curMZ;
        [ mVec ] = computeMapping( possibleRelatedPeaks', inpeakMZ, errRange );
        curPeak = mVec( ~isnan( mVec(:, 2) ), 2 );
        if isempty( DTemplate )
            DTemplate = zeros( SLEN, 1 );
            DTemplate(curPeak, 1) = 1;
            DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
            speciesM{cnt} = (possibleRelatedPeaks(1)-mpAry(1));
        else
            inFlag = 0;
            curTem = zeros( SLEN, 1 );
            curTem( curPeak, 1 ) = 1;
            for k = size( DTemplate, 2 ):-1:1
                tmp = DTemplate(:, k) - curTem;
                ins = sum(abs(tmp));
                if ins ==0
                    inFlag = 1;
                    speciesM{k} = [speciesM{k} possibleRelatedPeaks(1)-mpAry(1)];
                    DIonName{k} = [DIonName{k}  nameAry( ~isnan( mVec(:,2) ) )];
                end
            end
            if inFlag == 0
                cnt = cnt + 1;
                DTemplate = [DTemplate curTem];
                DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
                speciesM{cnt} = (possibleRelatedPeaks(1)-mpAry(1));
            end
        end
    end
end

if removeRedFlag == 1
    fprintf( 'computing redundancy in the DTemplate...\n' );
    %redundent is defined as:
    % A:[1 1 0 0]
    % B:[1 1 1 0]
    % A is redudent of B.
    redunVec = zeros( size( DTemplate, 2 ), 1 );
    % recVec = zeros( size(DTemplate, 2), 1);
    nonZeroEnt = cell( size(DTemplate, 2), 1 );
    nonZeroLen = zeros( size(DTemplate, 2), 1 );
    for i = 1:size(DTemplate, 2)
        nonZeroEnt{i} = find( DTemplate(:, i) == 1 );
        nonZeroLen(i) = length(nonZeroEnt{i});
    end
    
    for i = 1:size(DTemplate, 2)
        fprintf('%d\n', i );
        if i == 47
            keyboard
        end
        cNZLen = nonZeroLen(i);
        cEle = DTemplate(:, i);
        targetEle = find( nonZeroLen > cNZLen );
        for j = 1:length(targetEle)
            oEle = DTemplate(:, targetEle(j));
            ins = cEle | oEle;
            if length(find(ins==1)) == nonZeroLen(targetEle(j))
                redunVec(i) = 1;
                break;
            end
        end
    end
    tmp1 = DIonName;
    tmp2 = speciesM;
    DIonName = {};
    speciesM = {};
    cnt = 1;
    for i = 1:length(tmp1)
        if redunVec(i) == 0
            DIonName{cnt} = tmp1{i};
            speciesM{cnt} = tmp2{i};
            cnt = cnt + 1;
        end
    end
    DTemplate = DTemplate(:, redunVec==0);
end
end

