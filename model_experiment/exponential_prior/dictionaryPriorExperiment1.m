function dictionaryPriorExperiment1( )
%% First part: data generation for MASS dictionary, with exponential dropping
%data generation parameters
factor = 2;
lowerVal = 0.01;
highVal = 1;
dataNum = 1e6;
resolution = 1e4;

for i = 1:20
    if 0.5 / (factor^i) < lowerVal
        break;
    end
end
range = i; %plus it self
sampleMat = zeros( dataNum, range );

for j = 1:dataNum
    dominantPeak = rand(1,1)*0.5 + 0.5;
    peakList = zeros(range, 1);
    peakList(1) = dominantPeak;
    cPeak = dominantPeak;
    for i = 1:range-1
        cPeak = cPeak / factor;
        peakList(1+i) = cPeak;
    end
    sampleMat(j, :) = peakList;
end
figure;hist(sampleMat(:), resolution );

%% Second part: fit dist. to data
%matlab tool to fit data
dfittool;
% using bounded Pareto dist. to fit
logL = log(lowerVal);
sumLogX = sum(log(sampleMat(:)));
myfun = @(alpha)boundedParetoFunc( alpha, sampleMat(:), lowerVal, highVal, logL, sumLogX );
lb = zeros( 1, 1 );
ub = [];
options.MaxFunEvals = 1e4;
options.MaxIter = 1e4;
options.GradObj = 'on';
[lalpha,~] = fmincon( myfun,1e-3, [], [], [], [], lb, ub, [], options);
[nelements, centers]=hist(sampleMat(:), resolution);
totalN = sum(nelements)* (highVal-lowerVal) / resolution;
nelements = nelements / totalN;
figure;plot(centers,nelements);
hold on; plot( lowerVal:1e-4:highVal,boundedparetopdf( lowerVal, highVal, lalpha, lowerVal:1e-4:highVal ), 'r' );
legend('histogram of generated data','fitted Pareto distribution');
title( 'Fit Bounded-Pareto distribution to generated data, with alpha = 4.79e-9' );
end