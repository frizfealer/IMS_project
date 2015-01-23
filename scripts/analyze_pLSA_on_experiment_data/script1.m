%load input variables for negative mode data.
load( 'D:\Users\YeuChern\GitHub\IMS_project\example\realExp_100831_348_136_12_40hr_0-1XLB_LN\100831_348_136_12_40hr_0_1XLB_LN_input.mat' );

[~, pIdx1] = min(abs(mzAxis-1007.7));
[~, pIdx2] = min(abs(mzAxis-1021.8));
[~, pIdx3] = min(abs(mzAxis-1035.8));
% A possible error peak
[~, pIdx4] = min(abs(mzAxis-1057.6));

%positive peaks
[~, targetIdx(1)] = min(abs(mzAxis-1031));
[~, targetIdx(2)] = min(abs(mzAxis-1045));
[~, targetIdx(3)] = min(abs(mzAxis-1059));


%% shown in the data
% figure; subplot( 2, 2, 1); imagesc( reshape( dataCube(pIdx1, :), 38, 60 ) ); colorbar ;title('m/z = 1007.7');
% subplot( 2, 2, 2 ); imagesc( reshape( dataCube(pIdx2, :), 38, 60 ) ); colorbar ;title('m/z = 1007.7');
% subplot( 2, 2, 3 ); imagesc( reshape( dataCube(pIdx3, :), 38, 60 ) ); colorbar ;title('m/z = 1007.7');
% subplot( 2, 2, 4 ); imagesc( reshape( dataCube(pIdx4, :), 38, 60 ) ); colorbar ;title('m/z = 1007.7');

%% get target compounds in dictionary 
targetM = [];
for i = 1:size( nDTemplate, 2)
    if nDTemplate(pIdx1, i) == 1
        targetM = [targetM, i];
    elseif nDTemplate(pIdx2, i) == 1
        targetM = [targetM, i];
    elseif nDTemplate(pIdx3, i) == 1
        targetM = [targetM, i];
    elseif nDTemplate(pIdx4, i) == 1
        targetM = [targetM, i];
    end
end
initD = nDTemplate(:, targetM);

targetM2 = [];
for i = 1:size( nDTemplate, 2)
    if nDTemplate(pIdx1, i) == 1 && nDTemplate(pIdx2, i) == 1 && nDTemplate(pIdx3, i) == 1
        targetM2 = [targetM2, i];
    end
end
initD2 = nDTemplate(:, targetM2);


%% test on pLSA
inY = dataCube(:, :);
D_pLSA_Cell = {};
LL_Cell = {};
AICVec = zeros( size(inY, 1), 1 );

Learn.Verbosity = 0;
Learn.Max_Iterations = 200; 
Learn.heldout = 0; % for tempered EM only, percentage of held out data
Learn.Min_Likelihood_Change = 1e-3;
Learn.Folding_Iterations = 20; % for TEM only: number of fiolding
                              % in iterations
Learn.TEM =0; %tempered o not tempered
% start pLSA

for i = 1:90
    %Par = []; Par.maxit = 1e3; Par.Leps = 1; Par.doplot = 0;
    %[Pw_z,Pd_z,Pz,Li] = pLSA_EM( inY(:, BlkDS.indMap==1), i, Par );
    [Pw_z,Pz_d,Pd,Li,perp,beta] = pLSA( sparse(inY(:, BlkDS.indMap==1)),[],i,Learn);
    for j = 1:i
        Pw_z(:, j) = Pw_z(:, j) / norm( Pw_z(:, j) );
    end
    D_pLSA_Cell{i} = Pw_z;
    LL_Cell{i} = Li(end);
    [ AICVec(i) ] = computeAICc( inY, Li(end), Pw_z, BlkDS, i );
    save( 'D_pLSA_LL_AIC_LP.mat', 'D_pLSA_Cell', 'LL_Cell', 'AICVec' );
    % figure; subplot( 2, 2, 1); plot( Pw_z(pIdx1, :) ); title('m/z = 1007.7');
    % subplot( 2, 2, 2 ); plot( Pw_z(pIdx2, :) ); title('m/z = 1021.8');
    % subplot( 2, 2, 3 ); plot( Pw_z(pIdx3, :) ); title('m/z = 1035.8');
    % subplot( 2, 2, 4 ); plot( Pw_z(pIdx4, :) ); title('m/z = 1057.6');
end

%% show LL for each K
LLVec = zeros( length( LL_Cell ), 1 );
for i = 1:length(LL_Cell)
    LLVec(i) = LL_Cell{i};
end
figure; plot(LLVec); xlabel( 'compound number' ); ylabel( 'log-likelihood' );
export_fig k_ll.pdf -transparent
%% show AICc for each K
figure; plot(AICVec(1:length(LLVec))); xlabel( 'compound number' ); ylabel( 'AICc' );
export_fig k_AICc.pdf -transparent
%% check peaks
statVec = zeros( length(LL_Cell), 1 );
for i = 1:150
    [ resVec ] = evaluateRecovery( D_pLSA_Cell{i}, [pIdx1, pIdx2, pIdx3], setdiff(1:396, [pIdx1, pIdx2, pIdx3]), 1e-2 );
    statVec(i) = max(resVec);
end
statVec = zeros( length(LL_Cell), 1 );
for i = 1:length(LL_Cell)
    [ resVec ] = evaluateRecovery( D_pLSA_Cell{i}, [pIdx1, pIdx2, pIdx3], setdiff(1:317, [pIdx1, pIdx2, pIdx3]), 1e-2 );
    statVec(i) = max(resVec);
end
    [ resVec ] = evaluateRecovery( D_pLSA_Cell{106}, [pIdx1, pIdx2, pIdx3], setdiff(1:317, [pIdx1, pIdx2, pIdx3]), 1e-2 );

%positive peaks
statVec = zeros( length(LL_Cell), 1 );
for i = 1:length(LL_Cell)
    [ resVec ] = evaluateRecovery( D_pLSA_Cell{i}, targetIdx, setdiff(1:355, targetIdx), 1e-2 );
    statVec(i) = max(resVec);
end
[ resVec ] = evaluateRecovery( D_pLSA_Cell{90}, targetIdx, setdiff(1:355, targetIdx), 1e-2);
[ ~, maxE ] = max( resVec );
figure; plot( D_pLSA_Cell{90}(:, maxE) ); xlabel( 'm/z value' ); ylabel( 'intensity' );
for i = 10:10:90
    figure; subplot( 2, 2, 1); plot( D_pLSA_Cell{i}(pIdx1, :) ); title('m/z = 1007.7');
    subplot( 2, 2, 2 ); plot( D_pLSA_Cell{i}(pIdx2, :) ); title('m/z = 1021.8');
    subplot( 2, 2, 3 ); plot( D_pLSA_Cell{i}(pIdx3, :) ); title('m/z = 1035.8');
    subplot( 2, 2, 4 ); plot( D_pLSA_Cell{i}(pIdx4, :) ); title('m/z = 1057.6');
end

figure; subplot( 2, 2, 1); plot( D_pLSA_Cell{end}(targetIdx(1), :) ); title('m/z = 1007.7');
subplot( 2, 2, 2 ); plot( D_pLSA_Cell{end}(targetIdx(2), :) ); title('m/z = 1021.8');
subplot( 2, 2, 3 ); plot( D_pLSA_Cell{end}(targetIdx(3), :) ); title('m/z = 1035.8');
subplot( 2, 2, 4 ); plot( D_pLSA_Cell{end}(pIdx4, :) ); title('m/z = 1057.6');

