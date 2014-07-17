function dictionaryPriorExperiment2()
%% parameters
sampleSize = 1;
specSize = 100;
elementSize = 50;

%% data without any prior
W = abs( randn( elementSize, sampleSize ) );
D = abs( randn( specSize, elementSize ) );
Y = D * W + abs( randn( specSize, sampleSize ) );
Y2 = D * W;
%non-noise learning
learnedW = zeros( size( W ) );
for i = 1:sampleSize
    learnedW(:,i) = lsqnonneg(D, Y2(:,i)); 
end
norm( learnedW - W )
learnedD = zeros( size( D ) );
Y2r = Y2';
for i = 1:specSize
    cY = Y2r(:,i);
    learnedD(i, :) = lsqnonneg(W',cY)';
end
norm( learnedD-D )
%noise learning
learnedW = zeros( size( W ) );
for i = 1:sampleSize
    learnedW(:,i) = lsqnonneg(D, Y(:,i)); 
end
norm( learnedW - W )
learnedD = zeros( size( D ) );
Yr = Y';
for i = 1:specSize
    cY = Yr(:,i);
    learnedD(i, :) = lsqnonneg(W',cY)';
end
norm( learnedD-D )
%% data with exponential prior
factor = 1.4;
lowerVal = 0.01;
nFold = 5;

cVal = 0.8;
for i = 1:20
    fprintf( '%g\n', cVal / factor^i);
    if 0.5 / (factor^i) < lowerVal
        break;
    end
end
maxRan = i + 1;
W = zeros( elementSize, sampleSize );
for i = 1:sampleSize
    a = 0.5 + (0.5).*rand(1,1);
    ran = randi( maxRan-4, 1, 1 ) + 4;
    idx = sort( randi(elementSize,ran,1) );
    W(idx,i) = a*1.414.^-([0:1:(ran-1)]);
%     W(idx,i) = abs( randn(length(idx), 1) );
end
Y = D * W + abs( randn( specSize, sampleSize ) ); 
learnedW = zeros( size( W) );
for i = 1:sampleSize
    learnedW(:,i) = lsqnonneg( D,Y(:, i) );
end
norm( Y - D*learnedW ).^2/size(Y,1)

% lambdaGrid = 0:0.001:0.3;
% CVErrGridAry = zeros( length( lambdaGrid ), 1 );
% trainErrGridAry = zeros( length( lambdaGrid ), 1 );
% indices = crossvalind('Kfold', specSize, nFold );
opts=[];
opts.init=2;        % starting from a zero point
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations
opts.nFlag=0;       % without normalization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
opts.rsL2=0;        % the squared two norm term
% for z = 1:length(lambdaGrid)
%     learnedW2 = zeros( size( W ) );
%     CVErrAry = zeros( nFold, 1 );
%     trainErrAry = zeros( nFold, 1 );
%     for i = 1:nFold
%         test = (indices == i); train = ~test;
% %         param.lambda = lambdaGrid(z);
% %         param.lambda2 = 0;
% %         param.numThreads=-1;
% %         param.mode=2;
% %         param.pos = true;
%         for j = 1:sampleSize
% %             [~, path]=mexLasso(cY,cW,param);
%             [x, ~]= nnLeastR(D(train,:), Y(train,j), lambdaGrid(z), opts);
%             learnedW2(:,j) = x;
%         end
%         residual = Y(test, :) - D(test, :)*learnedW2;
%         trainErr = Y(train,:)-D(train,:)*learnedW2;
%         trainErrAry(i) = sum(trainErr(:).^2)/length(train);
%         CVErrAry(i) = sum(residual(:).^2)/length(test);
%     end
%     trainErr = sum(trainErrAry)/nFold;
%     trainErrGridAry(z) = trainErr;
%     CVErr = sum(CVErrAry)/nFold;
%     CVErrGridAry(z) = CVErr;
% end
%  [~,idx]=min(CVErrGridAry);
% param.lambda = lambdaGrid(idx);
% param.lambda2 = 0;
% param.numThreads=-1;
% param.mode=2;
% param.pos = true;
learnedW2 = zeros( size( W ) );
for j = 1:sampleSize
%     [~, path]=mexLasso(cY,cW,param);
%     learnedD2(j,:) = path(:,end)';
%     [B,FitInfo] = lasso( cW, cY, 'CV', 5 );
     [x, ~]= nnLeastR(D, Y(:, j), 0.2, opts);
    learnedW2(:, j) =x;
end
% figure; subplot(3,1, 1); plot(W);   
% subplot(3,1, 2); plot(learnedW);
% subplot(3,1,3); plot(learnedW2);

learnedW3 = ones( size( W ) )*1e-3;
alpha =  1e-3;
beta = 0.01;
L = 0.001;
H = 1;
logAlpha = log( alpha );
logLAlpha = log(L)*alpha;
thirdT = log(1-(L/H)^alpha);
options.GradObj = 'on';
for i = 1:sampleSize
myfun = @(x)exponetialPriorFunc( x, Y(:, i), D, logAlpha, logLAlpha, thirdT, alpha, beta );
x0 = learnedW3(:, i);
lb = zeros( elementSize, 1 );
ub = ones( elementSize, 1 );
[x,fval] = fmincon( myfun,x0,[],[],[],[], lb, ub, [], options);
learnedW3(:, i) = x;
end

figure; subplot(4,1, 1); plot(W); title( 'ground truth W' );
subplot(4,1, 2); plot(learnedW); title( 'learned W without prior' );
subplot(4,1,3); plot(learnedW2); title( 'W with L1 prior, lambda = 0.2' );
subplot(4, 1, 4); plot(learnedW3); title ('W with bounded Pareto prior, alpha = 1e-3, beta = 0.01' );

end