function [ outD ] = initD( dataCube, DTemplate, method )
%initD initialize dictioanry
mLen = size( DTemplate, 2 );
sLen = size( DTemplate, 1 );
% thres = round( mLen*0.05 );
thres = 10;
outD = DTemplate;
if strcmp( method, 'random' ) == 1
    for i = 1:mLen
        cSpec =  DTemplate(:, i) ~= 0 ;
        ins = dataCube(cSpec, :);
        res = sum(ins, 1);
        [~,idx]= sort( res, 'descend' );
        if find(res(idx)<=0,1) < thres
            thres = find(res(idx)<=0,1) - 1;
        end
        idx = idx(1:thres);
        dIdx = idx(randi(thres));
        outD(cSpec, i) = ins(:, dIdx);
        if norm( outD(cSpec, i), 2 ) == 0
            keyboard;
            dIdx = idx(randi(thres));
            outD(cSpec, i) = ins(:, dIdx);
        end
        outD(:, i) = outD(:, i)./ norm( outD(cSpec, i), 2 );
    end
elseif strcmp( method, 'online-DL' ) == 1
    fprintf( 'under construction...\n' );
%     y = dataCube(:, :);
%     y = log(y);
%     y(y ==-inf) = 0;
%     param.mode = 2; %min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2 + ...
%     % lambda||alpha_i||_1 + lambda_2||alpha_i||_2^2
%     param.modeD = 0; %C={ D in Real^{m x p} s.t. forall j, ||d_j||_2^2 <= 1 }
%     param.D = DTemplate;
%     param.lambda = 0; param.lambda2 = 0;
%     param.posAlpha = true;
%     param.posD = true;
%     param.numThreads=-1; % number of threads
%     param.iter=50;
%     param.gamma1=0.3
%     D = mexTrainDL( y,param );
elseif strcmp( method, 'NNMF' ) == 1
     for i = 1:mLen
        cSpec =  DTemplate(:, i) ~= 0 ;
        ins = dataCube(cSpec, :);
        res = sum(ins, 1);
        [~,idx]= sort( res, 'descend' );
        lowerIdx = find(res(idx)<=0,1) - 1;
        [w,h]=nnmf( ins(:,idx(1:lowerIdx)), 1 );
        outD(cSpec, i) = w;
        outD(:, i) = outD(:, i)./ norm( outD(cSpec, i), 2 );
    end
end


end

