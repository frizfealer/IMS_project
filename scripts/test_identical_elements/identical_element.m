%% construct data
D = zeros(2, 2);
D(1,:) = 1;
sNum = 400;
W = abs(randn(2, sNum));
W0 = zeros( sNum, 1);
preY = D*W;
Y = poissrnd( exp(preY) );
%% setting parameters
lambda = 0.5; theta = 0;
z0 = log(Y);
z0(z0==-inf)=0;
BlkDS.indMap = ones(400, 1);
aMatrix = ones(400, 1);
%if not in traing set, we shoud not initialize z0 according to it.
%take the average values
tmp = sum( z0(:, :), 2 ) / ( length( find( BlkDS.indMap == 1 ) ) );
for i = 1:sNum
    %if is in data area, but not in training set, 
    %set it to average spectrum
    if BlkDS.indMap(i) && ~aMatrix(i)
        z0(:, i) = tmp;
    end
end
scaleFactor =  1 / ( max( Y(:) ) / max( z0(:) ) ) * 10;

W_init = zeros( 2, sNum );
W0_init = zeros( sNum, 1 );
z1 = W_init;
BlkDS.blkNum = 1;
BlkDS.B2GMap{1} = 1:sNum;
[Rall] = genSparseGroupingMatrix( sNum, 1, 1 );
for i = 1:BlkDS.blkNum
    Rblk{i} = [];
    Rb = Rall(:, BlkDS.B2GMap{i});
    ins = sum(abs(Rb), 2);
    Rb(ins<2,:) = [];
    Rblk{i} = Rb;
end
%% first test, with zero initialization on W and same value on D
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
export_fig -transparent same_D_zero_W_init.pdf
close
%% second test, with one element random initialization on W and same value on D
W_init(1,:) = abs(randn(1,sNum));
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
export_fig -transparent same_D_zero_W_random_init.pdf
close
%% third test, with zero initialization on W and different value on D
D(1, 1) = 0.9; D(1, 2) = 0.3;
W_init = zeros( 2, sNum );
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
WCell0 = WResStruct.WHistCell;
export_fig -transparent same_D_diff_W_zero.pdf
close
%% forth test, with zero initialization on W and little different value on D
D(1, 1) = 0.55; D(1, 2) = 0.5;
W_init = zeros( 2, sNum );
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
WCell1 = WResStruct.WHistCell;
export_fig -transparent same_D_diff_0.55_0.5_W_zero.pdf
close
%% fivth test, with zero initialization on W and little different on D and lambda = 0.75
lambda = 0.75;
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
WCell2 = WResStruct.WHistCell;
export_fig -transparent same_D_diff_0.55_0.5_W_zero_lambda_0.75.pdf
close
%% sixth test, with zero initialization on W and litte different on D and lambda = 0.5
D(1, 1) = 0.9; D(1, 2) = 0.85; lambda = 0.5;
W_init = zeros( 2, sNum );
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
WCell3 = WResStruct.WHistCell;
export_fig -transparent same_D_diff_0.9_0.85_W_zero_lambda_0.5.pdf
close
%% seventh test, with zero initialzation on W and large different on D lambda = 0.1
D(1, 1) = 0.9; D(1, 2) = 0.3;
W_init = zeros( 2, sNum );
lambda = 0.1;
[WResStruct] = updateW_ADMM( Y, D, W_init, W0_init, z0, z1, aMatrix, BlkDS, 100, lambda, 0, theta, scaleFactor, [], Rblk, 1 );
figure; subplot(1,2,1); plot(WResStruct.W(1,:)); title('first element''s weight'); subplot(1,2,2); plot(WResStruct.W(2,:)); title('second element''s weight' );
WCell4 = WResStruct.WHistCell;
export_fig -transparent same_D_diff_0.9_0.3_W_zero_lambda_0.1.pdf
close

%% show the changes of W among iterations
figure;
for i = 1:100
    imagesc(WCell2{i});title(['i = ', num2str(i)]); colorbar;
    pause(0.75)
end

