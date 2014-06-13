function [ rSet ] = estimateDRelevant( dataCube, outD, percent )
%estimateRelevant estimate relevant of each dictionary elements
if isempty( percent )
    percent = 0.2;
end
mLen = size( outD, 2 );
Y = dataCube(:, :);
Y = log(Y); Y(Y==-inf)=0;
resAry = zeros( mLen, 1 );
parfor i = 1:mLen
    wi = Y \ outD(:, i);
    res = Y - outD(:, i)*wi';
    res = res(:);
    resAry(i) = sum( res.^2 );
end
[~,idx] = sort(resAry);
rSet = idx( 1:round( mLen*percent ) );
end

