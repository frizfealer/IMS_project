function [ gD, uDTemplate ] = synthesizeD( type, params, verbose )
%synthesizeD synthesize dictionary
%type: 'DTemplate'
%params: DTemplate, 
%   sparsePrec: the percentage of entries used in each
%   dictionary element. 
%   coheMax: the coherence max allowed of each element
%type: 'random'
%params: SLEN, MLEN
%   sparsePrec: the percentage of entries in all entries. e.g. SLEN*MLEN
%   coheMax: the coherence max allowed of each element
    
if strcmp( type, 'DTemplate' )
    assert( isfield( params, 'DTemplate' ) );
    assert( isfield( params, 'sparsePrec' ) );
    assert( isfield( params, 'coheMax' ) );
    sparsePrec = params.sparsePrec; coheMax = params.coheMax; DTemplate = params.DTemplate;
    [SLEN, MLEN] = size( DTemplate );
    %% generate current dictionary template with sparseVal condition
    redIdx = [];
    uDTemplate = zeros( SLEN, MLEN );
    for i = 1:MLEN
        if verbose == 1
            fprintf( 'element %d\n', i );
        end
        %     moleculePAry = round( appearingIonMS(i)+mpAry )'; %M+H, M+Na, ...etc.
        %     moleculePAry(moleculePAry<1) = [];
        %     moleculePAry(moleculePAry>sLen) = [];
        %     DTemplate(moleculePAry,i) = 1;
        moleculePAry = find( DTemplate(:,i) > 0 );
        usedNum = floor( length( moleculePAry ) * sparsePrec );
        if usedNum == 0
            usedNum = 1;
        end
        usedMZ = sort( randperm( length( moleculePAry ), usedNum ) );
        uDTemplate(moleculePAry(usedMZ), i) = 1;
        for j = 1:i-1
            if norm( uDTemplate(:,i)-uDTemplate(:,j), 1 ) == 0
                redIdx = [redIdx, i];
                break;
            end
        end
    end
    redIdx = unique( redIdx );
    uDTemplate(:,redIdx) = [];
    rMLEN = size( uDTemplate, 2 );
    gD = zeros( SLEN, rMLEN );
    
    numVec = zeros( rMLEN, 1 );
    for i = 1:rMLEN
        numVec(i) = length( find( uDTemplate(:, i) == 1 ) );
    end
    [~, numIdx]=sort( numVec );
    %% construct m/z map to species
    specMap = cell( SLEN, 1 );
    for i = 1:SLEN
        specMap{i} = find( uDTemplate(i, :)~=0 );
    end
    %% construct current dictionary template with coheMax condition
    availIdx = ones( rMLEN, 1 );
    for i = numIdx'
        if verbose == 1
            fprintf( 'element: %d\n', i );
        end
        infElements = [];
        moleculePAry = find( uDTemplate(:,i) > 0 );
        for j = moleculePAry'
            tmp = specMap{j};
            tmp = tmp( availIdx(tmp) == 0 );
            infElements = [infElements tmp];
        end
        infElements = unique( infElements );
        availIdx(i) = 0;
        gD(moleculePAry, i) = abs( randn( length( moleculePAry ), 1 ) );
        gD(:,i) = gD(:,i)/norm( gD(:,i), 2);
        if ~isempty( infElements )
            for j = infElements
                alpha = 1e-1;
                %fprintf( 'cov with %d: %g\n', j, gD(:,i)'*gD(:,j) );
                while abs( gD(:,i)'*gD(:,j) ) > coheMax
                    if abs( gD(:,i)'*gD(:,j) ) > 1
                        keyboard();
                    end
                    %fprintf( 'change of cov: %g\n', gD(:,i)'*gD(:,j) );
                    moleculePAry2 = find( uDTemplate(:,j) > 0 );
                    targetPAry = intersect( moleculePAry, moleculePAry2 );
                    gD(targetPAry, i) = gD(targetPAry, i) - alpha;
                    gD(gD(:, i)<0, i) = 0;
                    alpha = alpha + 0.01;
                end
            end
        end
        %         temp = gD(:,i)'*gD(:,1:i-1);
        %         %fprintf( 'debug: %g\n', max(temp) );
        %         if max(temp) > coheMax
        %             keyboard();
        %         end
    end
elseif strcmp( type, 'random' ) == 1
    assert( isfield( params, 'SLEN' ) );
    assert( isfield( params, 'MLEN' ) );
    assert( isfield( params, 'sparsePrec' ) );
    assert( isfield( params, 'coheMax' ) );
    sparsePrec = params.sparsePrec; coheMax = params.coheMax;
    SLEN = params.SLEN; MLEN = params.MLEN;
%% generate current dictionary template with sparseVal condition
    usedNum = floor( SLEN*MLEN * sparsePrec );
    usedEnt = sort( randperm( SLEN*MLEN, usedNum ) );
    uDTemplate = zeros( SLEN, MLEN );
    uDTemplate(usedEnt) = 1;
    for i = 1:MLEN
        if norm( uDTemplate(:, i) )== 0
            uDTemplate( randi(MLEN), i ) = 1;
        end
    end
    numVec = zeros( MLEN, 1 );
    for i = 1:MLEN
        numVec(i) = length( find( uDTemplate(:, i) == 1 ) );
    end
    [~, numIdx]=sort( numVec );
    %% construct m/z map to species
    specMap = cell( SLEN, 1 );
    for i = 1:SLEN
        specMap{i} = find( uDTemplate(i, :)~=0 );
    end
    %% construct current dictionary template with coheMax condition
    gD = zeros( SLEN, MLEN );
    availIdx = ones( MLEN, 1 );
    for i = numIdx'
        if verbose == 1
            fprintf( 'element: %d\n', i );
        end
        infElements = [];
        moleculePAry = find( uDTemplate(:, i) > 0 );
        for j = moleculePAry'
            tmp = specMap{j};
            tmp = tmp( availIdx(tmp) == 0 );
            infElements = [infElements tmp];
        end
        infElements = unique( infElements );
        availIdx(i) = 0;
        gD(moleculePAry, i) = abs( randn( length( moleculePAry ), 1 ) );
        gD(:,i) = gD(:,i)/norm( gD(:,i) );
        if ~isempty( infElements )
            for j = infElements
                alpha = 1e-1;
                %fprintf( 'cov with %d: %g\n', j, gD(:,i)'*gD(:,j) );
                while abs( gD(:,i)'*gD(:,j) ) > coheMax
                    if abs( gD(:,i)'*gD(:,j) ) > 1
                        keyboard();
                    end
                    %fprintf( 'change of cov: %g\n', gD(:,i)'*gD(:,j) );
                    moleculePAry2 = find( uDTemplate(:,j) > 0 );
                    targetPAry = intersect( moleculePAry, moleculePAry2 );
                    if norm(gD(:, i), 2 ) > norm( gD(:, j), 2 )
                        gD(targetPAry, i) = gD(targetPAry, i) - alpha;
                    else
                        gD(targetPAry, j) = gD(targetPAry, j) - alpha;
                    end
                    gD(gD(:, i)<0, i) = 0;
                    alpha = alpha + 0.01;
                end
            end
        end
    end
end
    
end

