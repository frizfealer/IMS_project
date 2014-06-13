function [ gD, uDTemplate, usedTerm ] = synthesizeD( SLEN, MLEN, DTemplate, sparseVal, coheVal, verbose )
%synthesizeD synthesize dictionary
%sparseVal: 0~1, 1 means no sparse, at least used one signal in an element
assert( SLEN == size( DTemplate, 1 ) );
assert( MLEN <= size( DTemplate, 2 ) );

usedTerm = sort( randperm( size( DTemplate, 2 ), MLEN ) );
redI = [];
uDTemplate = DTemplate( :, usedTerm );

for i = 1:MLEN
    if verbose == 1
        fprintf( 'element %d\n', i );
    end
    %     moleculePAry = round( appearingIonMS(i)+mpAry )'; %M+H, M+Na, ...etc.
    %     moleculePAry(moleculePAry<1) = [];
    %     moleculePAry(moleculePAry>sLen) = [];
    %     DTemplate(moleculePAry,i) = 1;
    moleculePAry = find( uDTemplate(:,i) > 0 );
    usedNum = floor( length( moleculePAry ) * sparseVal );
    if usedNum == 0
        usedNum = 1;
    end
    usedMZ = sort( randperm( length( moleculePAry ), usedNum ) );
    uDTemplate(:, i) =0;
    uDTemplate(moleculePAry(usedMZ), i) = 1;
    for j = 1:i-1
        if norm( uDTemplate(:,i)-uDTemplate(:,j), 1 ) == 0
            redI = [redI, i];
        end
    end
end
redI = unique( redI );
uDTemplate(:,redI) = [];


rMLEN = size( uDTemplate, 2 );
gD = zeros( SLEN, rMLEN );
specMap = cell( SLEN, 1 );
for i = 1:SLEN
    specMap{i} = find( uDTemplate(i, :)~=0 );
end

for i = 1:rMLEN
    if verbose == 1
        fprintf( 'element: %d\n', i );
    end
    infElements = [];
    moleculePAry = find( uDTemplate(:,i) > 0 );
    for j = moleculePAry'
        tmp = specMap{j};
        tmp = tmp( tmp < i );
        infElements = [infElements tmp];
        infElements = unique( infElements );
    end
    gD(moleculePAry, i) = abs( randn( length( moleculePAry ), 1 ) );
    gD(:,i) = gD(:,i)/norm( gD(:,i), 2);
    if ~isempty( infElements )
        for j = infElements
            alpha = 1e-1;
            %fprintf( 'cov with %d: %g\n', j, gD(:,i)'*gD(:,j) );
            while abs( gD(:,i)'*gD(:,j) ) > coheVal
                if abs( gD(:,i)'*gD(:,j) ) > 1
                    keyboard();
                end
                %fprintf( 'change of cov: %g\n', gD(:,i)'*gD(:,j) );
                moleculePAry2 = find( uDTemplate(:,j) > 0 );
                targetPAry = intersect( moleculePAry, moleculePAry2 );
                gD(targetPAry, i) = gD(targetPAry, i) - alpha;
                gD(gD(:, i)<0, i) = 0;
                alpha = alpha + 0.02;
            end
        end
    end
    temp = gD(:,i)'*gD(:,1:i-1);
    %fprintf( 'debug: %g\n', max(temp) );
    if max(temp) > coheVal
        keyboard();
    end
end

    
end

