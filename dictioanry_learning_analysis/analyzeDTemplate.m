function [ statStruc ] = analyzeDTemplate( DTemplate, statFlag )
%analyze_DTemplate analyze DTemplate (or dictionary) generatd from molecular patterns
%statFlag = [0,1], is a vector of 0, 1 that 0 means donot do the test and 1
%means does
[sLen, mLen] = size( DTemplate );
if statFlag(1) == 1
    statAry1 = zeros( sLen, 1 );
    for i = 1:sLen
        statAry1(i) = length( find( DTemplate(i,:) ~= 0 ) );
    end
    figure; hist( statAry1, 25 ); title( 'histogram of non-zero-entries of dictionary elements' );
    xlabel( '# non-zero-entries' ); ylabel( 'counts' );
    statStruc.statAry1 = statAry1;
end
if statFlag(2) == 1
    statAry2 = zeros( size(DTemplate, 2), 1 );
    for i = 1:size(DTemplate, 2)
        fprintf( 'element %d\n', i );
        cTemplate = DTemplate(:, i);
        for j = (i+1):size(DTemplate, 2)
            ccTemplate = DTemplate(:, j);
            res = norm( cTemplate - ccTemplate, 1 );
            if res == 0
                statAry2(i) = statAry2(i)+1;
            end
        end
    end
    figure; hist( statAry2, 25 ); 
    title( 'histogram of redundant dictioanry elements' );
    xlabel( 'redundant number' ); ylabel( 'counts' );
    statStruc.statAry2 = statAry2;
end

if statFlag(3) == 1
    statAry3 = zeros( mLen, 1 ); 
    for i = 1:mLen
        statAry3(i) = norm( DTemplate(:, i), 2 );
    end
    figure; hist( statAry3 ); title( 'histogram of norm of dictioanry elements' );
    xlabel( 'norm of dictioanry elements' ); ylabel( 'counts' );
    statStruc.statAry3 = statAry3;
end

if statFlag(4) == 1
    statAry4 = zeros( mLen, 1 );
    for i = 1:mLen
        temp = DTemplate(:,i)'*DTemplate;
        temp(i) = [];
        statAry4(i) = max( temp );
    end
    figure; hist( statAry4 ); title( 'histogram of corr. of dictioanry elements' );
    xlabel( 'correlations of dictionary elements' ); ylabel( 'counts' );    
    statStruc.statAry4 = statAry4;
end

end

