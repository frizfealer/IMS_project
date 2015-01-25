function [ DTemplate, DIonName, speciesM ] = genDTemplate( inpeakMZ, MPPATH, smallestM, errRange, varargin )
%genDTemplate generate DTemplate within epsilon error
%inpeakMZ [sLen 1],
%MPPATH: the pattern file path,
%'D:\RA\Project_dictionaryLearning_IMS\molecule_profile_new_trim.csv'
%smallestM: the smallest molecule's weight
%errRange: the error tolerance to bin m/z to inpeakMZ
%inMolMZ: the additional parameter that users can specified the m/z to
%generate peaks, instead of using the whole inpeakMZ
%return DTemplate, a template with 1 means possible and 0 means impossible
%DIonName, "possible" adducts generated that signals. This information is
%just one of the possible answers.
fileID = fopen(MPPATH);
C = textscan(fileID, '%f %s', 'delimiter', {','});
fclose(fileID);
mpAry = C{1,1};
nameAry = C{1,2};
SLEN = length(inpeakMZ);
if length(varargin) == 1
    inMolMZ = varargin{1};
else
    inMolMZ = inpeakMZ;
end
% DMSIndex = [];

% interval = mpAry - mpAry(1);
% mapping = zeros( length(interval), length(interval) );
% for i = 2:length(interval)
%     mapping(i,:) = interval - interval(i);
% end
% mapping(1,:) = interval;

MVec = [];
fprintf( 'generating possible MWs...\n' );
for i = 1:length(inMolMZ)
%     fprintf( 'm/z %g\n', inMolMZ(i));
        for j = 1:length(mpAry)
            if inMolMZ(i) >= (smallestM+mpAry(j))
                MVec = [MVec inMolMZ(i)-mpAry(j)];
            end
        end
end
MVec = sort(MVec);
curPos = 1;
mMVec = MVec(curPos);
for i = 2:length(MVec)
    if MVec(i) - MVec(curPos) > errRange
        curPos = i;
        mMVec = [mMVec MVec(i)];
    end  
end

speciesM = [];
DTemplate = [];
DIonName = {};
cnt = 1;
fprintf( 'generating DTemplate...\n' );
for i = 1:length(mMVec)
%     fprintf( 'm/z %g\n', mMVec(i));
    possibleRelatedPeaks = mpAry+mMVec(i);
    [ mVec ] = computeMapping( possibleRelatedPeaks, inpeakMZ, errRange );
    curPeak = mVec( ~isnan( mVec(:, 2) ), 2 );
    if isempty( DTemplate )
        DTemplate = zeros( SLEN, 1 );
        DTemplate(curPeak, 1) = 1;
        DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
        speciesM = [speciesM mMVec(i)];
    else
        curTem = zeros( SLEN, 1 );
        curTem( curPeak, 1 ) = 1;
        inFlag = 0;
        tmp = mMVec(i)-mpAry(end);
        %get the lower index in DTemplate we should compare for the newest
        %one
        lowerM = mMVec(find( mMVec >= tmp, 1 ));
        lowerIdx = find( speciesM >= lowerM );
        for k = size( DTemplate, 2 ):-1:lowerIdx
            tmp = DTemplate(:, k) - curTem;
            ins = sum(abs(tmp));
            if ins ==0
                inFlag = 1;
                break;
            end
        end
        if inFlag == 0
            cnt = cnt + 1;
            DTemplate = [DTemplate curTem];
            DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
            speciesM = [speciesM mMVec(i)];
        end
    end
end

fprintf( 'computing redundancy in the DTemplate...\n' );
%redundent is defined as:
% A:[1 1 0 0]
% B:[1 1 1 0]
% A is redudent of B.
redunVec = zeros( size( DTemplate, 2 ), 1 );
% recVec = zeros( size(DTemplate, 2), 1);
for i = 1:size(DTemplate, 2)
    cSet = DTemplate(:,i);
    for j = i:size(DTemplate, 2)
        tmp = cSet - DTemplate(:,j);
        if max(tmp)-min(tmp)==1 && min(tmp) == -1
            redunVec(i) = 1;
%             recVec(i) = j;
            break;
        end
    end
end
tmp = DIonName;
DIonName = {};
cnt = 1;
for i = 1:length(tmp)
    if redunVec(i) == 0
        DIonName{cnt} = tmp{i};
        cnt = cnt + 1;
    end
end
DTemplate = DTemplate(:, redunVec==0);
speciesM = speciesM(redunVec == 0);
% for i = 1:length(inMolMZ)
%     fprintf( 'm/z %g\n', inMolMZ(i));
%     for j = 1:length(mpAry)
%         if inMolMZ(i) >= (smallestM+mpAry(j))
%             possibleRelatedPeaks = mapping(j,:)+inMolMZ(i);
%             [ mVec ] = computeMapping( possibleRelatedPeaks', inpeakMZ', errRange );
%             curPeak = mVec( ~isnan( mVec(:, 2) ), 2 );
%             if isempty( DTemplate )
%                 DTemplate = zeros( SLEN, 1 );
%                 DTemplate(curPeak, 1) = 1;
%                 DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
%                 speciesM = [speciesM inMolMZ(i)-mpAry(j)];
%             else
%                 curTem = zeros( SLEN, 1 );
%                 curTem( curPeak, 1 ) = 1;
%                 inFlag = 0;
%                 for k = size( DTemplate, 2 ):-1:1
%                     ins = norm( DTemplate(:, k) - curTem, 1 );
%                     if ins ==0
%                         inFlag =1;
%                         break;
%                     end
%                 end
%                 if inFlag == 0
%                     cnt = cnt + 1;
%                     DTemplate = [DTemplate curTem];
%                     DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
%                     speciesM = [speciesM inMolMZ(i)-mpAry(j)];
%                 end
%                 %                     ins = DTemplate - repmat( curTem, 1, size( DTemplate, 2 ) );
%                 %                     ins = sum( abs(ins), 1);
%                 %                     if isempty( find( ins==0, 1 ) )
%                 %                         cnt = cnt + 1;
%                 %                         DTemplate = [DTemplate curTem];
%                 %                         DIonName{cnt} = nameAry( ~isnan( mVec(:,2) ) );
%                 %                     end
%             end
%         else
%             break;
%         end
%     end
% end


end
%
% testing
% for i = 1:length(DIonName)
%     if size(DIonName{i}, 1) ~= length(find(DTemplate(:,i)==1))
%         fprintf( '%d\n', i );
%     end
% end
