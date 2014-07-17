mLen = size( expRec.outW, 1 );
corrAry = zeros( mLen, 1 );
%grp1 and grp2
for i = 1:mLen
    t1 = full( expRec.outW(i,BlkDS.B2GMap{1}) );
    t2 = full( expRec.outW(i,BlkDS.B2GMap{3}) );
    wVal = max(max(t1), max(t2));
    fprintf( '%g\n', max(max(t1), max(t2)));
    if wVal < 1
        corrAry(i) = 0;
    else
        h = ttest2(t1,t2);
        corrAry(i) = h;
    end
end
idx = find( corrAry==1 );
ins = idx(end-30:end);
for i = 1:31
    figure; imagesc(reshape(expRec.outW(ins(i),:),54,88));colorbar;title(['index= ', num2str(ins(i))]);
end