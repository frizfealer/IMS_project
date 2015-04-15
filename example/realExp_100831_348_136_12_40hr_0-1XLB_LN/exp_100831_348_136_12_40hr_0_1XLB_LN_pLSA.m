function [ output_args ] = exp_100831_348_136_12_40hr_0_1XLB_LN_pLSA( inVarPath, outResPath )
% inVarPath = '/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/100831_348_136_12_40hr_0_1XLB_LN_vars.mat';
% outResPath = '/csbiohome01/ycharn/IMS_project-master/example/realExp_100831_348_136_12_40hr_0-1XLB_LN/res_pLSA_20150403.mat';
load(inVarPath);
testNum = 50;
Pw_zCell = cell(testNum, 1); AICVec = zeros( testNum, 1 ); AICVec2 = zeros( testNum, 1 );
LiVec = zeros( testNum, 1);
parfor i = 1:testNum
    [~, Pw_zCell{i}, ~, ~, Li] = pLSA_EM_simple( IMSD.dataCube(:, IMSD.BlkDS.indMap==1)+1e-8, [], i );
    samNum = IMSD.dataCube(:, IMSD.BlkDS.indMap==1); samNum = length(samNum(:));
    idx = find(Li, 1, 'last');
    [AICVec(i)] = computeAICc_sparse(samNum, Li(idx), Pw_zCell{i}); 
    AICVec2(i) = computeAICc( IMSD.dataCube, LiVec(i), IMSD.BlkDS, i );
    LiVec(i) = Li(idx);
end
save( outResPath, 'Pw_zCell', 'AICVec', 'AICVec2', 'LiVec' );
%% run sparse-pLSA
load( outResPath );
load( inVarPath );
figure;plot(LiVec);
%% chose the component number to be the least number that has likelihood saturated.
bestIdx = 35;
[Pz_dw, ~, ~, ~, ~] = pLSA_EM_simple( IMSD.dataCube(:, IMSD.BlkDS.indMap==1)+1e-8, [], bestIdx );
ssVec = 0:0.00000005:0.0000005;
AICSPVec = zeros( length(ssVec), 1 );
Pw_zSPVec = cell(length(ssVec), 1); Pd_zSPVec = cell( length(ssVec), 1 );
parfor i = 1:length(ssVec)
    [Pw_zSPVec{i}, Pd_zSPVec{i}, ~, ~, LL] = pLSA_EM_sparse( IMSD.dataCube(:, IMSD.BlkDS.indMap==1)+1e-8, bestIdx, 0, ssVec(i), Pz_dw );
    samNum = length(find(IMSD.BlkDS.indMap==1));
    [AICSPVec(i)] = computeAICc_sparse(samNum, LL, Pw_zSPVec{i} );
end
save( outResPath, 'Pw_zCell', 'AICVec', 'AICVec2', 'LiVec', 'Pw_zSPVec', 'Pd_zSPVec', 'AICSPVec' );
end

