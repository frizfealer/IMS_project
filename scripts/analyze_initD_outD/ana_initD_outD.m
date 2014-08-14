%A script to analyze the difference between Dinit (D comes from
%initialization) and outD( D comes from dictionary learning with the same
%initialization)
%% A example of generate subplot automatically
MLEN = size( pDTemplate, 2 );
ins = length( MLEN, 1 );
for i = 1:MLEN
    ins(i) = norm( expRec.outD(:, i) );
end
target = find( ins <= 0.1 );
FIG_W_NUM = 4;
FIG_H_NUM = 4;
prevWinNum = 1;
figure;
for i = 1:length(target)
    cWinNum = floor( (i-1) / (FIG_W_NUM*FIG_H_NUM) ) + 1;
    if cWinNum ~= prevWinNum
        prevWinNum = cWinNum;
        figure;
    end
    cyPos = floor( (i-1) / FIG_W_NUM )+1 - FIG_H_NUM * (cWinNum - 1);
    cxPos = floor( mod( (i-1), FIG_W_NUM ) ) + 1;
    fprintf('%d, %d %d\n',cWinNum, cyPos, cxPos );
    cSpec = find(pDTemplate(:,target(i))==1);
    subplot( FIG_H_NUM, FIG_W_NUM, (cyPos-1)*FIG_W_NUM+cxPos ); 
    plot( 1:length(cSpec), DFirst(cSpec, target(i)), 'b-o', 1:length(cSpec), expRec.outD(cSpec, target(i)), 'g-*' );
    legend('initial D', 'out D'); 
    title( ['L2-norm = ', num2str( norm( DFirst(:, target(i))-expRec.outD(:,target(i)) ) )] );
end