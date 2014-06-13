function [ D ] = genSparseGroupingMatrix( hei, wid, repNum )
%genSparseGroupingMatrix generate sparse grouping matrix for 2D fused lasso
%term, used by updatez1_m
%D, [wid*(hei-1)+hei*(wid-1) hei*wid]
% %build up E matrix, the relationship between the Wm in the grid
% %The |w_i,j-w_i-1,j| relationship
% repNum: the repeat number, the number of different molecular species z1
% in my cases
interval = hei*wid;
yAry1 = zeros( 2, wid*(hei-1) );
xAry1 = zeros( 2, wid*(hei-1) );
dataAry1 = ones( 2, wid*(hei-1) );
dataAry1(1,:) = -1;
cnt = 1;
for i=1:(interval-1)
    if mod( i, hei ) ~= 0
        yAry1(:, cnt) = i;
        xAry1(1, cnt) = i;
        xAry1(2, cnt) = i + 1;
        cnt = cnt + 1;
    end
end
yAry1 = yAry1(:);
xAry1 = xAry1(:);
dataAry1 = dataAry1(:);
E1 = sparse( yAry1, xAry1, dataAry1, interval, interval );
% nrow = [];
% for i = 1:size( E1, 1 )
%     if ismember( 1, E1( i, : ) )
%         nrow = [nrow; i];
%     end
% end
% E1 = E1(nrow, :);
%The |w_i,j-w_i,j-1| relationship
interval = hei*wid;
yAry2 = zeros( 2, hei*(wid-1) );
xAry2 = zeros( 2, hei*(wid-1) );
dataAry2 = ones( 2, hei*(wid-1) );
dataAry2(1,:) = -1;
for i = 1:(interval-hei)
    yAry2(:,i) = i;
    xAry2(1,i) = i;
    xAry2(2,i) = i+hei;
end
yAry2 = yAry2(:);
xAry2 = xAry2(:);
dataAry2 = dataAry2(:);
E2 = sparse( yAry2, xAry2, dataAry2, interval, interval );
% nrow = [];
% for i = 1:size( E2, 1 )
%     if ismember( 1, E2( i, : ) )
%         nrow = [nrow; i];
%     end
% end
% E2 = E2(nrow, :);
% D = [E1; E2];
cYAry = yAry2 + interval;
cYAry = [yAry1; cYAry];
cXAry = [xAry1; xAry2];
cDataAry = [dataAry1; dataAry2];
% D = sparse( cYAry, cXAry, cDataAry, interval*2, interval );

len = length( cYAry );
aYAry = zeros( len*repNum, 1 );
aXAry = zeros( len*repNum, 1 ); 
aDataAry = zeros( len*repNum, 1);
for i = 1:repNum
    aYAry( ((i-1)*len+1):(i*len) ) = cYAry+(i-1)*interval*2;
    aXAry( ((i-1)*len+1):(i*len) ) = cXAry+(i-1)*interval;
    aDataAry( ((i-1)*len+1):(i*len) ) = cDataAry;
end
D = sparse( aYAry, aXAry, aDataAry, interval*2*repNum, interval*repNum );

%remove all-zero lines
ins = sum(abs(D),2);
res = ins==0;
D(res,:)=[];
end
