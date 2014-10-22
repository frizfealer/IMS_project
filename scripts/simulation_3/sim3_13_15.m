function [ output_args ] = sim3_13_15( )
load( 'simulation3_env_1.mat' );
addpath( genpath( '~/IMS_project-master' ) );
addpath( genpath( '~/SLEP_package_4.1/SLEP' ) );
addpath( '~/Ipopt_3118' );
param.CLUSTER_NAME = 'killdevil1024_1'; %usually this is the default name
param.CLUSTER_NUM = 12; %number of cluster uses, usually 12
resVec = zeros(1:27, 2);
expTemp = [];
for i = 13:15
    lambda = fVec(i);
    theta = sVec(i);
    phi = tVec(i);
    vvv = zeros(5,1);
    for j = 1:5
        if j == 1
            param.OUTER_IT_NUM = 30;
            initVal = [];
        else
            param.OUTER_IT_NUM = 20;
            initVal = expRec;
        end
        aMatrix = ones(30, 30);
        aMatrix(foldMap(:)==j)=0;
        [ expRec ] = dictionaryLearning_ADMM_v5( pY, initVal, sDTemplate, [], lambda, theta, phi, aMatrix, pBlkDS, [], snapPath, param );
        expTemp{i,j} = expRec;
        save( 'resTemp_2.mat', 'expTemp', '-v7.3' );
        ttt=aMatrix==0;
        [ vvv(j) ] = validationOnTesting( ttt, pY, expRec.outW, expRec.outW0, expRec.outD );
    end
    resVec(i, 1) = mean(vvv);
    resVec(i, 2) = std(vvv);
    save( 'resVecTemp_2.mat', 'resVec' );
end
end


