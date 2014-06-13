function [ z1Idx, COMNUM ] = genSquareBlockIdx( CUNIT, hei, wid )
%genSquareBlockIdx generate block indexes
%for dividing the input grid into squares
%input:
%CUNIT (Computing Unit) the length of the square
%hei, wid, the size of the input grid
%output:
%z1Idx, cell with length = COMNUM, number of the CUNIT
%for each unit i z1Idx{i} is a [m 3] matrix, where m = number of nodes in
%that units
% [m, 1] [m, 2] is the [height, width] coordinates of that nodes,
% [m, 3] is the transformed coordinates of that nodes, (2D into 1D, by
% column), this is done by height+(width-1)*height
cnt = 1; %area counter
% CUNIT = 10;
COMNUM = ( floor( (hei - 1 ) / CUNIT ) + 1 )*( floor( ( wid - 1) / CUNIT ) + 1 );
z1Idx = cell( COMNUM, 1);
for i = 1:CUNIT:hei
    for j = 1:CUNIT:wid
        if i == 1
            startH = 1;
        else
            startH = i - 1;
        end
        if i + CUNIT - 1 > hei
            endH = hei;
        else
            endH = i + CUNIT - 1;
        end
        if j == 1
            startW = 1;
        else
            startW = j - 1;
        end
        if j + CUNIT -1 > wid
            endW = wid;
        else
            endW = j + CUNIT - 1;
        end
        cnt2 = 1;
        z1Idx{cnt} = zeros( (endH-startH+1)*(endW-startW+1), 3);
        for b = startW:endW
            for a = startH:endH
                z1Idx{cnt}(cnt2,1:2) = [a b];
                z1Idx{cnt}(cnt2, 3) = a + (b-1)*hei;
                cnt2 = cnt2 + 1;
            end
        end
        cnt = cnt + 1;
    end
end

outMap = zeros(hei, wid);
for j = 1:COMNUM
    for i = 1:size(z1Idx{j}, 1)
%         outMap( z1Idx{j}(i, 1),z1Idx{j}(i, 2) ) = outMap( z1Idx{j}(i, 1),z1Idx{j}(i, 2) ) + j;
        outMap( z1Idx{j}(i, 3) ) = outMap( z1Idx{j}(i, 3) ) + j;
    end
end
figure;imagesc(outMap);
end

