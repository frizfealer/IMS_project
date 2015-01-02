function [ scaleFactor ] = computeScaleFactor( Y, aMatrix )
%computeScaleFactor 
%compute scaleFactor
BlkDS = conBLKDS(Y);
[~, hei, wid] = size(Y);
tmp = log(Y); tmp(tmp==-inf)=0; 
idx = intersect( find(aMatrix==1), find(BlkDS.indMap==1) );
tmp2 = sum( tmp(:, idx), 2 ) / ( length( idx ) );
for i = 1:hei*wid
    if BlkDS.indMap(i) == 1 && aMatrix(i) == 0
        tmp(:, i) = tmp2;
    end
end
scaleFactor =  1 / ( max( Y(:) ) / max( tmp(:) ) );

end

