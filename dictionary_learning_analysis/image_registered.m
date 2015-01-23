function [fHandle, mapBioSam, MMapMSSam, macSam ] = image_registered( bioFilePath, macFilePath, MSFilePath, requireInfoPath, alphaVal, verbose )
% bioFilePath='registration_example\DSC01573.jpg';
% macFilePath='registration_example\DSC00000.jpg';
% MSFilePath='registration_example\example.jpg';
% requireInfoPath='registration_example\registerInfo.mat'
% shiftFlag=1, for this example. Because the DSC00000 is differet from
% DSC01573
% alphaVal, the value of transparent, 0 means non-visible.
% verbose = 1, generating intermediate results
bioSam = imread( bioFilePath );
MSSam = imread( MSFilePath );
macSam = imread( macFilePath );

load( requireInfoPath );
mapBioSam = imwarp( bioSam, mytForm2 );
mapMSSam = imwarp( MSSam, mytForm );
[hei, wid, ~] = size( mapMSSam );
MMapMSSam = zeros( hei, wid, 3, 'uint8' );
MMapMSSam(1:(hei-shiftHEI), 1:(wid-shiftWID), :) = mapMSSam((shiftHEI+1):hei, (shiftWID+1):wid, :);
if verbose == 1
    figure;imshowpair( macSam, mapBioSam );
    figure;imshowpair( macSam, MMapMSSam );
end
%     figure;imshowpair( rMSSam, rrBioSam );
% figure('Visible','Off'); 
imshow( mapBioSam ); hold on; h = imshow( MMapMSSam ); 
% indMap = zeros( size( rMSSam, 1 ), size( rMSSam, 2 ) );
% for i = 1:size( rMSSam, 1  )
%     for j = 1:size( rMSSam, 2 )
%         if rMSSam(i, j, 1) ~= 0 || rMSSam(i, j, 2) ~= 0 || (rMSSam(i, j, 3) ~= 142 )
%             indMap(i, j) = 1;
%         end
%     end
% end
% alphaMat = zeros( size( rMSSam, 1 ), size( rMSSam, 2 ) );
% alphaMat(indMap==1) = alphaVal;
alphaMat = zeros( size(MMapMSSam, 1), size(MMapMSSam, 2) );
%20 is a fixed threshold, depending on the colormap of MS image.
idx1 = ( MMapMSSam(:,:,1)>=20 | MMapMSSam(:,:,2) >= 20 | MMapMSSam(:,:,3) >= 20 );
alphaMat(idx1)=alphaVal;
set( h, 'AlphaData', alphaMat );
fHandle = gcf;
end
