function model_selection_deal_Data( dataCube, expRecReal, pDTemplate, pSpeciesM, pDIonName, mzAxis )
%this function is deprecated, seems the degree of freedom is hard to
%estimate in fusion terms
BlkDS = conBLKDS( dataCube );
[ sigDict, sigM, sigIdx] = analyzeResults( expRecReal, pDTemplate, BlkDS, pSpeciesM, pDIonName, mzAxis, [], 892, '~/', 0 );
D = expRecReal.outD;
IGNORE_INT = 1e-2;
Y = sparse( dataCube(:, BlkDS.indMap==1) );

%% option of ipopt section      
% Set up the auxiliary data.

% Set the IPOPT options.
options.ipopt.print_level = 0;
% options.ipopt.jac_d_constant   = 'yes';
% options.ipopt.hessian_constant = 'yes';
options.ipopt.mu_strategy      = 'adaptive';
options.ipopt.max_iter         = 200;
options.ipopt.tol              = 1e-6;
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.limited_memory_max_history = 1500; 
options.ipopt.honor_original_bounds = 'yes';
options.ipopt.bound_relax_factor = 0;

% The callback functions.
funcs.objective         = @objective;
funcs.gradient          = @gradient;
funcs.hessian          = @hessian;
funcs.hessianstructure = @hessianstructure;

%% test for two hypothesis for each dictionary elements
%one is one dictionary elements
%another is n dictionary elements (assume a dictionary element with n
%entries > IGNORE_INT), nd each elements only has one entry.
for i = 1:size(sigDict, 2)
    res(i) = struct( 'val1', 0, 'val2', 0, 'val_overDis1', 0, 'val_overDis2', 0, 'c_hat1', 0, 'c_hat2', 0, 'pval', 0 );
end

LLAry = zeros( 2, size(D, 2 ) );
for i = 1:size(D, 2)
    cMZ = find( D(:, sigIdx(i)) > IGNORE_INT );
    firstD = D(cMZ, sigIdx(i));
    firstD = sparse( [firstD ones(length(cMZ), 1)] );
    startWF = zeros( 2, 1 );
    firstW = zeros( 2, size(Y, 2) );
    % The constraint functions are bounded from below by zero.
    options.lb = zeros( length(startWF),1 );
    options.ub = Inf( length(startWF), 1 );
    tic
    parfor j = 1:size(Y, 2)
        o{j} = options;
        o{j}.auxdata = { Y(cMZ, j), firstD };
        [x, ~] = ipopt_auxdata(startWF,funcs,o{j});
        firstW(:, j) = x;
    end
    toc
    
    secondD = zeros( length(cMZ) );
    for j = 1:length(cMZ)
        secondD(j, j) = D(cMZ(j), sigIdx(i));
    end
    secondD = sparse( [secondD ones(length(cMZ), 1)] );
    secondW = zeros( (length(cMZ)+1), size(Y, 2) );
    startWS = zeros( (length(cMZ)+1), 1 );
    options.lb = zeros( length(startWS),1 );
    options.ub = Inf( length(startWS), 1 );
    tic
    parfor j = 1:size(Y, 2)
        o{j} = options;
        o{j}.auxdata = { Y(cMZ, j), secondD };
        [x, ~] = ipopt_auxdata(startWS,funcs,o{j});
        secondW(:, j) = x;
    end
    toc
    LLAry(1, i) = Poisson_LL_Func( Y(cMZ, :), firstD, firstW, [] );
    LLAry(2, i) = Poisson_LL_Func( Y(cMZ, :), secondD, secondW, [] ); 
    save( 'res_LLAry.mat', 'LLAry' );
%     [ val1, val_overDis1, c_hat1 ] = compueAICC( @Poisson_LL_Func, Y(cMZ, :), firstD, firstW, length(firstW(:)) );
%     [ val2, val_overDis2, c_hat2 ] = compueAICC( @Poisson_LL_Func, Y(cMZ, :), secondD, secondW, length(secondW(:)) );
%     pval = computePVal_PoissonReg( Y(cMZ, :), firstD, firstW, secondD, secondW, length(firstW(:)), length(secondW(:)) );
%     res(i).val1 = val1; res(i).val2 = val2;
%     res(i).val_overDis1 = val_overDis1; res(i).val_overDis2 = val_overDis2;
%     res(i).c_hat1 = c_hat1; res(i).c_hat2 = c_hat2;
%     res(i).pval = pval;
end
save( 'res_temp.mat', 'res' );

%%use res_LLAry to test the result
Y = sparse( dataCube(:, BlkDS.indMap==1) );
D = expRecReal.outD;
load('res_LLAry.mat');
logY = log(Y);
logY(isinf(logY)) = 0;
LLY = Poisson_LL_Func( Y, [], [], logY );
pureAIC = zeros(2, 30 );
AICC = zeros(2, 30);
AICC_od = zeros(2, 30);
c_hat = zeros(2, 30);
pVal = zeros(1, 30 );
[ sigDict, sigM, sigIdx] = analyzeResults( expRecReal, pDTemplate, BlkDS, pSpeciesM, pDIonName, mzAxis, [], 892, '~/', 0 );
IGNORE_INT = 1e-2;
for i = 1:30
     cMZ = find( D(:, sigIdx(i)) > IGNORE_INT );
     p = 1*size(Y,2);
     q = length(cMZ)*size(Y,2);
     n = q;
    [ pureAIC(1,i), AICC(1,i), AICC_od(1, i), c_hat(1, i) ] = compueAICC( LLAry(1, i), LLY, p, n  );
    [ pureAIC(2,i), AICC(2,i), AICC_od(2, i), c_hat(2, i) ] = compueAICC( LLAry(2, i), LLY, q, n );
    pVal(i) = computePVal_PoissonReg( LLAry(1, i), LLAry(2, i), p, q );
end
out = zeros( 2*30, 5 );
for i = 1:30
    start = (i-1)*2;
    out((start+1):(start+2), 1) = pureAIC(:, i);
    out((start+1):(start+2), 2) = AICC(:, i);
    out((start+1):(start+2), 3) = AICC_od(:, i);
    out((start+1):(start+2), 4) = c_hat(:, i);
    out((start+1):(start+2), 5) = [pVal(i); 0];
end

end

function f = objective(x, auxdata)
%x is current W
[Y, firstD] = deal(auxdata{[1 2]});
samSize = size(Y, 2);
x = reshape(x, size(firstD, 2), samSize );
preY = firstD*x;
f = sum( sum( -Y.*preY + exp(preY) ) );
end

function g = gradient(x, auxdata)
%x is current W
[Y, firstD] = deal(auxdata{[1 2]});
samSize = size(Y, 2);
x = reshape(x, size(firstD, 2), samSize );
ePreY = exp(firstD*x);
g = -firstD'*Y + firstD'*ePreY;
end