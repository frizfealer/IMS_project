function [] = exp_100831_348_136_12_40hr_0_1XLB_LN_NNMF( inVarPath, outResPath )
% inVarPath = '/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_vars.mat';
% outResPath = '/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/res_NNMF_20150403.mat';
load(inVarPath);
D_NNMF_Cell = {};
ResVec = [];
for i = 1:50
    rNum = 1e2; mNum = i;
    [D_NNMF_Cell{i}, ResVec(i)] = NNMF_als_wrapper( IMSD.dataCube, [], IMSD.BlkDS, mNum, rNum, true, 'default' );
end
save( outResPath, 'D_NNMF_Cell', 'ResVec');
%run sparse-NNMF
load( outResPath);
load( inVarPath );
% [~, best_NNMF] = min( ResVec );
best_NNMF = 17;
option.eta = 0; option.beta = 0;
betaVec = 0:0.1:3;
ResSpVec = zeros( length(betaVec), 1 );
D_NNMF_SP_Cell = cell( length(betaVec), 1 );
option.iter = 500;
parfor i = 1:length(betaVec)
    opt = option;
    opt.beta = betaVec(i);
    [D_NNMF_SP_Cell{i},~,~,~,ResSpVec(i)] = sparsenmfnnls( IMSD.dataCube(:, IMSD.BlkDS.indMap==1 ), best_NNMF, opt );
end
save( outResPath, 'D_NNMF_Cell', 'ResVec', 'D_NNMF_SP_Cell', 'ResSpVec');

% betaVec2 = 7:15;
% ResSpVec2 = zeros( length(betaVec), 1 );
% D_NNMF_SP_Cell2 = cell( length(betaVec), 1 );
% option.iter = 500;
% parfor i = 1:length(betaVec2)
%     opt = option;
%     opt.beta = betaVec2(i);
%     [D_NNMF_SP_Cell2{i},~,~,~,ResSpVec2(i)] = sparsenmfnnls( IMSD.dataCube(:, IMSD.BlkDS.indMap==1 ), best_NNMF, opt );
% end


% save('NNMF_res_v2.mat', 'D_NNMF_Cell', 'ResVec', 'D_NNMF_SP_Cell', 'ResSpVec');

end

