function [ outD ] = initD( dataCube, DTemplate, method, DIonName, L1Flag )
%--------------------------------------------------------------------------
% initD: initialize dictionary according to DTemplate
%--------------------------------------------------------------------------
% DESCRIPTION:
%   initialize dictionary according to dictionary template, with different
%   algorithms. Available algortihms are 'random' and 'NNMF'.
%
% INPUT ARGUMENTS:
%   dataCube, size of [s w h], that is #(m/z)*width*height
%   DTemplate, size of [s m], a dictionary pattern.
%   method, 'random', 'online-DL' (under construction), 'NNMF'
%   DIonName, can be [], if provided, CH2-ion-type is set to zero,
%   deprecated.
%   L1Flag, should be always set to 0 (Under construction).
% OUTPUT ARGUMENTS:
%   outD, the resulting dictionary.

mLen = size( DTemplate, 2 );
% thres = round( mLen*0.05 );
thres = 10;
outD = zeros( size( DTemplate ) );
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
        cSpec =  find( DTemplate(:, i) ~= 0 );
        if ~isempty( DIonName )
            cIonName = DIonName{i};
            ch2Idx = [];
            for j = 1:length(cIonName)
                if ~isempty( strfind( cIonName{j}, 'CH2' ) )
                    ch2Idx = [ch2Idx j];
                end
            end
            cSpec(ch2Idx) = [];
        end
        ins = dataCube(cSpec, :);
        res = sum(ins, 1);
        [~,idx]= sort( res, 'descend' );
        lowerIdx = find(res(idx)<=0,1) - 1;
        if isempty( lowerIdx )
            lowerIdx = length(idx);
        end
        if lowerIdx ~= 0
            [w, ~] = nnmf( ins(:,idx(1:lowerIdx)), 1 );
            %                 if isempty( DIonName )
            %                     outD(cSpec, i) = w;
            %                     outD(:, i) = outD(:, i)./ norm( outD(cSpec, i), 2 );
            %                 else
            %                     cIonName = DIonName{i};
            %                     for j = 1:length(cIonName)
            %                         if ~isempty( strfind( cIonName{j}, 'CH2' ) )
            %                             w(j) = 1e-3;
            %                         end
            %                     end
            %                     outD(cSpec, i) = w;
            %                     outD(:, i) = outD(:, i)./ max( norm( outD(cSpec, i), 2 ), 1 );
            %                 end
            outD(cSpec, i) = w;
            %                 outD(cSpec, i) = 1;
            if L1Flag == 0
                outD(:, i) = outD(:, i)./ norm( outD(cSpec, i), 2 );
            elseif L1Flag == 1
                outD(:, i) = outD(:, i)./ norm( outD(cSpec, i), 1 );
            end
        end
        
        %% display and debugging...
        %             [lambda,psi,T,stats,F] = factoran( ins(:,idx(1:lowerIdx))', 1);
        % %             lambda = lambda / norm(lambda)
        % %             w = w / norm(w)
        % %             nnn = mean( ins(:,idx(1:lowerIdx)), 2 );
        % %             nnn = nnn / norm(nnn)
        %             corrcoef( ins(:,idx(1:lowerIdx))' )
        %             temp = ins(:,idx(1:lowerIdx))';
        %             [sa, sb] = size( temp );
        %             [ temp, centers ] = discretize( temp(:), 1e5 );
        %             temp = reshape( temp, sa, sb );
        %             rrr = find( cSpec == 1 );
        %             val = zeros( length(rrr), length(rrr) );
        %             for j = 1:length(rrr)-1
        %                 for z = j+1:length(rrr)
        %                     [ val(j, z) ] = computeMI( temp(:,j), temp(:,z), centers );
        %                 end
        %             end
        %             val
        %             figure;
        %             for j = 1:length(rrr)
        %                 subplot(3,3,j); imagesc(reshape(dataCube(rrr(j),:),54,88)); colorbar;
        %             end
    end
end

end

% % testing if CH2 is not initialze > 1e-3
% for i = 1:size(pDTemplate, 2)
%     cIonName = pDIonName{i};
%     idx = find(pDTemplate(:,i)==1);
%     temp =[];
%     for j = 1:length(cIonName)
%         if ~isempty( strfind(cIonName{j}, 'CH2' ))
%             temp = [temp j];
%         end
%     end
%     ttt = find(D_init(idx, i) > 1e-3);
%     if ~isempty(intersect(ttt,temp))
%         fprintf( 'element %d\n', i );
%     end
% end
