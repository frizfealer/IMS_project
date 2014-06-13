peakMZ = h5read('~/Documents/MATLAB/RA/2013_Bmyc_Paeni_Early_LP.h5','/entry_0/analysis_0/peak_mz');
Y = h5read('~/Documents/MATLAB/RA/2013_Bmyc_Paeni_Early_LP.h5','/entry_0/analysis_0/peak_cube');
% peakMZ = h5read('D:\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP.h5','/entry_0/analysis_0/peak_mz');
% Y = h5read('D:\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP.h5','/entry_0/analysis_0/peak_cube');
[sLen, HEI, WID] = size(Y);
%set smallest ion to be H+, 1.007276
%analysis peak cube, remove < 1.007276+1.007276
idx = find( peakMZ >= 1.007276+1.007276 , 1 );
rY = zeros( size(Y, 1)-idx+1, size(Y, 2), size(Y, 3) );
for i = 1:size(Y, 2)*size(Y, 3)
    rY(:,i) = Y(idx:end, i);
end
rpeakMZ = peakMZ(idx:end);

% %analysis of the difference between m/z
% rpeakMZdiff = [rpeakMZ(2:end); 0]-rpeakMZ;
% hist(rpeakMZdiff(1:end-1));
% ins = find(rpeakMZdiff(1:end-1)<1 );
% hist(rpeakMZdiff(ins));title('difference of m/z that is smaller then 1'); 
% ins = find(rpeakMZdiff(1:end-1)>1 );
% hist(rpeakMZdiff(ins));title('difference of m/z that is larger then 1'); 

%% binning to integer values of orignal m/z +/- 0.5 Da
brpeakMZ = unique(round(rpeakMZ));
binY = zeros(length(brpeakMZ), HEI, WID); 
for i = 1:length(brpeakMZ)
    for j = 1:length(rpeakMZ)
        if rpeakMZ(j) <= brpeakMZ(i) + 0.5
            if rpeakMZ(j) >= brpeakMZ(i) - 0.5
                binY(i,:) = binY(i,:) + rY(j,:);
            end
        else
            break;
        end
    end
end

%% generate template
% MPPath = 'D:\RA\Project_dictionaryLearning_IMS\molecule_profile_new.csv';
MPPath = '/csbiohome01/ycharn/Documents/molecule_profile_new.csv';
mpAry = csvread( MPPath );
[ DTemplate, DMSIndex ] = genDTemplate( 0.5, brpeakMZ, mpAry, 1, 1:brpeakMZ(end) );
% for i = 1:1432
% if find( sum(DTemplate(:,i))==0 )
%     fprintf('%d\n', i );
% end
% end

%% get a portion of data, and test on this part of data
maxV = max(binY(:));
[x,y]=find(binY(:,:)==maxV);
%binY(168,6,39) higest value

insBY = binY(:,6:11,39:44);
aMatrix = ones(6,6);
initVar = [];
[ LPAry, outD, outW, outW0, u0, u1, u2, u3, u4, z0, z1, z2, z3, z4 ] = ...
            dictionaryLearning_ADMM( insBY, initVar, DTemplate, 32, 0, 32, [], aMatrix, 100 );
cW = outW;
 mCW = cW;
mCW(mCW<=1e-3) = 0;
[mLen, HEI, WID] = size(mCW);
nonzeroAry = zeros(HEI*WID,1);
for i = 1:HEI*WID
    nonzeroAry(i) = length(find(mCW(:,i)>0));
end

L1Mat = zeros(5,5);
for i = 1:6
    for j = 1:6
        L1Mat(i,j) = norm(cW(:,i)-cW(:,j),1);
    end
end

diffMat = zeros(6,6)
for i = 1:6
    for j = 1:6
        diffMat(i,j) = log(max( abs( insBY(:,i)-insBY(:,j))));
    end
end
        
        
[cSLen, cHEI, cWID ] = size( insBY );
% [ aMatrix ] = genHeldOutAlpha( cHEI, cWID, 0.1 );
save( '01032014_ins_aMatrix.mat', aMatrix);
initVar = [];
seq = 1:20;
seq1 = 2.^seq;
lambda = 1e-3*seq1;
theta = 1e-3*seq1;





%% run on 5:10
for i = 5:10
    for j = 1:length(theta)
        fprintf('%d %d\n', i, j);
        [ LPAry, outD, outW, outW0, u0, u1, u2, u3, u4, z0, z1, z2, z3, z4 ] = ...
            dictionaryLearning_ADMM( insBY, initVar, DTemplate, lambda(i), 1, theta(j), [], aMatrix, 100 );
        if i==5 && j == 1
            expResults_5_10.LPAry = LPAry; expResults_5_10.outD = outD; expResults_5_10.outW = outW;
            expResults_5_10.outW0 = outW0; expResults_5_10.u0 = u0; expResults_5_10.u1 = u1;
            expResults_5_10.u2 = u2; expResults_5_10.u3 = u3; expResults_5_10.u4 = u4;
            expResults_5_10.z0 = z0; expResults_5_10.z1 = z1; expResults_5_10.z2 = z2;
            expResults_5_10.z3 = z3; expResults_5_10.z4 = z4;
        else
            expResults_5_10(i,j).LPAry = LPAry; expResults_5_10(i,j).outD = outD; expResults_5_10(i,j).outW = outW;
            expResults_5_10(i,j).outW0 = outW0; expResults_5_10(i,j).u0 = u0; expResults_5_10(i,j).u1 = u1;
            expResults_5_10(i,j).u2 = u2; expResults_5_10(i,j).u3 = u3; expResults_5_10(i,j).u4 = u4;
            expResults_5_10(i,j).z0 = z0; expResults_5_10(i,j).z1 = z1; expResults_5_10(i,j).z2 = z2;
            expResults_5_10(i,j).z3 = z3; expResults_5_10(i,j).z4 = z4;
        end
        save( '01042014_insResult_lambda_5to10.mat', 'expResults_5_10');
    end
end

%% run on 11:15
for i = 11:15
    for j = 1:length(theta)
        fprintf('%d %d\n', i, j);
        [ LPAry, outD, outW, outW0, u0, u1, u2, u3, u4, z0, z1, z2, z3, z4 ] = ...
            dictionaryLearning_ADMM( insBY, initVar, DTemplate, lambda(i), 1, theta(j), [], aMatrix, 100 );
        if i==11 && j == 1
            expResults_11_15.LPAry = LPAry; expResults_11_15.outD = outD; expResults_11_15.outW = outW;
            expResults_11_15.outW0 = outW0; expResults_11_15.u0 = u0; expResults_11_15.u1 = u1;
            expResults_11_15.u2 = u2; expResults_11_15.u3 = u3; expResults_11_15.u4 = u4;
            expResults_11_15.z0 = z0; expResults_11_15.z1 = z1; expResults_11_15.z2 = z2;
            expResults_11_15.z3 = z3; expResults_11_15.z4 = z4;
        else
            expResults_11_15(i,j).LPAry = LPAry; expResults_11_15(i,j).outD = outD; expResults_11_15(i,j).outW = outW;
            expResults_11_15(i,j).outW0 = outW0; expResults_11_15(i,j).u0 = u0; expResults_11_15(i,j).u1 = u1;
            expResults_11_15(i,j).u2 = u2; expResults_11_15(i,j).u3 = u3; expResults_11_15(i,j).u4 = u4;
            expResults_11_15(i,j).z0 = z0; expResults_11_15(i,j).z1 = z1; expResults_11_15(i,j).z2 = z2;
            expResults_11_15(i,j).z3 = z3; expResults_11_15(i,j).z4 = z4;
        end
        save( '01042014_insResult_lambda_11to15.mat', 'expResults_11_15');
    end
end

%% validation on testing data
logFInsBY = insBY;
for i = 1:cHEI*cWID
%     fprintf('%d\n', i);
    logFInsBY(:,i) = logFac( insBY(:,i) );
end
testValAry = zeros( size(expResults_6_10) ); 
for i = 1:size(expResults_6_10, 1)
    for j=1: size(expResults_6_10, 2)
        if ~isempty(expResults_6_10(i,j).LPAry)
            cW = expResults_6_10(i,j).outW; cW0 = expResults_6_10(i,j).outW0; cD = expResults_6_10(i,j).outD;
            [ testValAry(i,j) ] = validationOnTesting( 1-aMatrix, insBY, cW, cW0, cD, logFInsBY );
        else
            testValAry(i,j) = -inf;
        end
    end
end

figure; plot(testValAry(1,:)); title(['lambda= ', num2str(lambda(6)), ' theta=x-values']); xlabel('theta values');ylabel('log-likelihood on testing set');
figure; plot(testValAry(2,:)); title(['lambda= ', num2str(lambda(7)), ' theta=x-values']);xlabel('theta values');ylabel('log-likelihood on testing set');
figure; plot(testValAry(3,:)); title(['lambda= ', num2str(lambda(8)), ' theta=x-values']);xlabel('theta values');ylabel('log-likelihood on testing set');
figure; plot(testValAry(4,:));  title(['lambda= ', num2str(lambda(9)), ' theta=x-values']);xlabel('theta values');ylabel('log-likelihood on testing set');
figure; plot(testValAry(5,:));  title(['lambda= ', num2str(lambda(10)), ' theta=x-values']);xlabel('theta values');ylabel('log-likelihood on testing set');

%% the largest log-posterior is at lambda = 0.256, theta=1.024
i=8-5; j=10;
cW = expResults_6_10(i,j).outW; cW0 = expResults_6_10(i,j).outW0; cD = expResults_6_10(i,j).outD; cLPAry = expResults_6_10(i,j).LPAry;

%% analysis of the dictionary 

preY = cD*cW(:,:);
hist(preY(:));

maxCD = zeros(size(cD, 2), 1);
for i = 1:size(cD, 2)
    maxCD(i) = max(cD(:,i));
end
hist(maxCD);
zeroIdx = find(maxCD==0);
nonZeroIdx = setdiff(1:1432, zeroIdx);
imagesc(cD(:,nonZeroIdx)'*cD(:,nonZeroIdx));
mCW = cW;
mCW(mCW<=1e-3) = 0;
[mLen, HEI, WID] = size(mCW);
nonzeroAry = zeros(HEI*WID,1);
for i = 1:HEI*WID
    nonzeroAry(i) = length(find(mCW(:,i)>0));
end

mCD = cD;
mCD(mCD <= 1e-3 ) = 0;
nonzeroAry = zeros(1432,1);
for i = 1:1432
    nonzeroAry(i) = length(find(mCD(:,i)>0));
end