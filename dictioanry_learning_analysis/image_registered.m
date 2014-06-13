% bioSam = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\DSC01573.jpg');
% baseOfMac = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\DSC01582.jpg');
% msSam = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\example.jpg');
% [nRow, nCol, ~] = size( msSam );
% bioSam=imresize(bioSam,[nRow,nCol]);
% baseOfMac=imresize(baseOfMac,[nRow,nCol]);
% imwrite(bioSam,'bioSam.jpg');
% imwrite(baseOfMac,'baseOfMac.jpg');
% figure, imshow(bioSam)
% figure, imshow(baseOfMac)
% figure; imshow(msSam);
% cpselect(bioSam, orthophoto);
% cpselect(msSam, orthophoto);
% cpselect(msSam, bioSam);
% save( 'registeredPnt.mat', 'input_points', 'base_points','input_points1', 'base_points1');

bioSam = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\bioSam.jpg');
baseOfMac = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\baseOfMac.jpg');
msSam = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\example.jpg');
cpselect(bioSam, baseOfMac);
cpselect(msSam, baseOfMac);
save( 'registeredPnt.mat', 'input_points', 'base_points','input_points1', 'base_points1');

% rBioSam=imresize(rBioSam,[nRow,nCol]);
figure;
baseOfMac = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\DSC01582.jpg');
imshow(baseOfMac);

hold on
h = imshow(rBioSam);
set( h, 'AlphaData', .7 );

rMsSam = imwarp(msSam, mytform2);
%rMsSam=imresize(rMsSam,[nRow,nCol]);
figure;
h=imshow(baseOfMac);
hold on
h = imshow(rMsSam);
set(h,'AlphaData',.5);
hold off

figure;
h=imshow(rMsSam);
hold on
h = imshow(baseOfMac);
set(h,'AlphaData',.7);
hold off
alpha(0.1);

% load( 'registeredPnt.mat' );
% mytform = fitgeotrans(input_points, base_points, 'affine');
% save( 'transInfo.mat', 'mytform' );
% 
% registered = imwarp(bioSam, mytform);
% figure; subplot(1,2,1); imshow(baseOfMac);
% subplot(1,2,2); imshow(registered);

% figure; subplot(1,2,1); imshow(registered);
% subplot(1,2,2); imshow(bioSam);

msSam = imread('D:\Users\YeuChern\Dropbox\unc\CS\RA\Project_dictionaryLearning_IMS\data\2013_Bmyc_Paeni_Early_LP\example.jpg');
figure, imshow(msSam)
cpselect(msSam, orthophoto);
save( 'registeredPnt2.mat', 'input_points1', 'base_points1');
 figure;
 imshow(baseOfMac); 
 hold on 
 h = imshow(msSam); 
 hold off
 alpha(0.8);
keyboard();


registered2s = imwarp(msSam, mytform2);
figure; subplot(1,2,1); imshow(baseOfMac);
subplot(1,2,2); imshow(registered2s);

figure;
imshow(orthophoto);
hold on
h = imshow(registered);
hold off
alpha(0.8);


unregistered2 = imread('DSC01573.jpg');
unregistered2s = imresize(unregistered2, 0.1275);


cpselect(unregistered2s, unregistereds);
save( 'registeredPnt.mat', 'input_points', 'base_points', 'input_points2', 'base_points2' );

load( 'registeredPnt.mat' );
mytform2 = fitgeotrans(input_points2, base_points2, 'affine');
save( 'transInfo.mat', 'mytform', 'mytform2' );

registered2s = imwarp(unregistered2s, mytform2);

 figure;
 imshow(unregistereds); 
 hold on 
 h = imshow(registered2s); 
 hold off
 alpha(0.8);
 
 registered3 = imwarp(registered2s, mytform);
 figure;
imshow(orthophoto);
hold on
h = imshow(registered3);
hold off
alpha(0.8);