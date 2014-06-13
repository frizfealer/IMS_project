function [ gD, uDTemplate ] = synthesizeD_woTem( SLEN, MLEN, supNum, coheVal, verbose )
%synthesizeD synthesize dictionary

assert(SLEN == MLEN*supNum-MLEN+1);
gD = zeros(SLEN, MLEN);
uDTemplate = zeros(SLEN, MLEN);
SVal = 1;
SCALE = 2;
val = zeros( supNum, 1 );
cVal = SVal;
for i = supNum:-1:1
    val(i) = cVal;
    cVal = cVal * SCALE;
end
startPos = 1;
for i = 1:MLEN
    uDTemplate(startPos:(startPos+supNum-1), i) = 1;
    gD(startPos:(startPos+supNum-1),i) = val;
    gD(:, i) = gD(:, i) / norm( gD(:, i), 2 );
    %next dictionary element starts at with one position overlapping.
    startPos = startPos + supNum - 1;
end

specMap = cell( SLEN, 1 );
for i = 1:SLEN
    specMap{i} = find( uDTemplate(i, :)~=0 );
end

for i = 1:MLEN
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
                if norm(gD(:, i), 2 ) > norm( gD(:, j), 2 )
                    gD(targetPAry, i) = gD(targetPAry, i) - alpha;
                else
                    gD(targetPAry, j) = gD(targetPAry, j) - alpha;
                end
                gD(gD(:, i)<0, i) = 0;
                gD(gD(:, j)<0, j) = 0;
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

