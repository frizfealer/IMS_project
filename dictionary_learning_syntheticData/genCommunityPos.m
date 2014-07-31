function [com] = genCommunityPos( comNum, IHEIGHT, IWIDTH  )
for i = 1:comNum
    %compute community's width
    com(i).wid = round( IWIDTH / ( comNum + abs( randn( 1,1 )*0.5 ) ) );
    %compute community's height
    com(i).hei = round( IHEIGHT / ( 0.5*abs( randn( 1,1 )+1 ) ) );
    %compute community's starting height
    com(i).startHei = round( IHEIGHT / 10 + rand(1,1)*2 );
    if com(i).startHei + com(i).hei - 1 > IHEIGHT
        com(i).hei  = IHEIGHT - com(i).startHei + 1 -round( IHEIGHT / 10 );
    end
    %compute community's starting width
    if i == 1
        com(i).startWid = round( IHEIGHT / 10 );
    else
        com(i).startWid = com(i-1).startWid+com(i-1).wid + round( IHEIGHT / 10 );
    end
    if com(i).startWid + com(i).wid - 1 > IWIDTH
        com(i).wid  = IWIDTH - com(i).startWid + 1  -round( IHEIGHT / 10 );
    end
    %construct the width, height coordinates of boundary points (8 points)
    com(i).pntWid(1:2) = sort( com(i).startWid + randi( com(i).wid, 2, 1 ) - 1 );
    while com(i).pntWid(2)-com(i).pntWid(1) <= 4
        com(i).pntWid(1:2) = sort( com(i).startWid + randi( com(i).wid, 2, 1 ) - 1 );
    end
%     if com(i).pntWid(2)-com(i).pntWid(1) <= 4
%         com(i).pntWid(2) = com(i).pntWid(1) + 5;
%     end
    com(i).pntWid(3:4) = com(i).startWid + com(i).wid - 1;
    com(i).pntWid(5:6) = sort( com(i).startWid + randi( com(i).wid, 2, 1 ) - 1, 'descend' );
    while com(i).pntWid(5)-com(i).pntWid(6) <= 4
        com(i).pntWid(5:6) = sort( com(i).startWid + randi( com(i).wid, 2, 1 ) - 1, 'descend' );
    end
%     if com(i).pntWid(5)-com(i).pntWid(6) <= 4
%         com(i).pntWid(5) = com(i).pntWid(6) + 5;
%     end
    com(i).pntWid(7:8) = com(i).startWid;
    com(i).pntHei(1:2) = IHEIGHT + 1 - com(i).startHei;
    com(i).pntHei(3:4) = sort( IHEIGHT + 1 - com(i).startHei - randi( com(i).hei, 2, 1 ) + 1 );
    while com(i).pntHei(4)-com(i).pntHei(3) <= 4
        com(i).pntHei(3:4) = sort( IHEIGHT + 1 - com(i).startHei - randi( com(i).hei, 2, 1 ) + 1 );
    end
%     if com(i).pntHei(4)-com(i).pntHei(3) <= 4
%         com(i).pntHei(4) = com(i).pntHei(3) + 5;
%     end
    com(i).pntHei(5:6) = IHEIGHT + 1 - com(i).startHei - com(i).hei + 1;
    com(i).pntHei(7:8) = sort( IHEIGHT + 1 - com(i).startHei - randi( com(i).hei, 2, 1 ) + 1, 'descend' );
    while com(i).pntHei(7) - com(i).pntHei(8) <= 4
        com(i).pntHei(7:8) = sort( IHEIGHT + 1 - com(i).startHei - randi( com(i).hei, 2, 1 ) + 1, 'descend' );
    end
%     if com(i).pntHei(7)-com(i).pntHei(8) <= 4
%         com(i).pntHei(7) = com(i).pntHei(8) + 5;
%     end
    com(i).pntWid(9) = com(i).pntWid(1);
    com(i).pntHei(9) = com(i).pntHei(1);
    %compute the correct order for the convex hull of the 8 points
    K = convhull( com(i).pntWid, com(i).pntHei );
    com(i).pntHei = com(i).pntHei(K); com(i).pntWid = com(i).pntWid(K);
    %the 
    [X,Y] = meshgrid( 1:IWIDTH, 1:IHEIGHT ); 
    IN = inpolygon(X,Y,com(i).pntWid, com(i).pntHei);
    com(i).areaIdx = IN;
end
end