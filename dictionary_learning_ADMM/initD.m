function [ outD ] = initD( dataCube, DTemplate, method )
%initD initialize dictioanry
mLen = size( DTemplate, 2 );
thres = round( mLen*0.05 );
outD = DTemplate;
if strcmp( method, 'random' ) == 1
    for i = 1:mLen
        cSpec =  DTemplate(:, i) ~= 0 ;
        ins = dataCube(cSpec, :);
        res = sum(ins, 1);
        [~,idx]= sort( res, 'descend' );
        idx = idx(1:thres);
        dIdx = idx(randi(thres));
        outD(cSpec, i) = ins(:, dIdx);
        outD(:, i) = outD(:, i)./ norm( outD(cSpec, i), 2 );
    end
end


end

