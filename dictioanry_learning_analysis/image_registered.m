function [fHandle, rrBioSam, rMSSam, macSam ] = image_registered( bioFilePath, macFilePath, MSFilePath, requireInfoPath, shiftFlag, alphaVal, verbose )
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
if isempty( requireInfoPath )
    cpselect( bioSam, macSam );
    cpselect( MSSam, macSam );
%     save( 'registeredPnt.mat', 'input_points', 'base_points','input_points2', 'base_points2');
    mytForm = fitgeotrans( input_points, base_points, 'projective' );
    mytForm2 = fitgeotrans( input_points2, base_points2, 'projective' );
%     save( 'transInfo.mat', 'mytForm', 'mytForm2' );
    requireInfo.input_points = input_points;
    requireInfo.base_points = base_points;
    requireInfo.input_points2 = input_points2;
    requireInfo.base_points2 = base_points2;
    requireInfo.mytForm = mytForm;
    requireInfo.mytForm2 = mytForm2;
    save( 'registerInfo.mat', 'requireInfo' );
else
    load( requireInfoPath );
    rBioSam = imwarp( bioSam, requireInfo.mytForm );
    rMSSam = imwarp( MSSam, requireInfo.mytForm2 );
end
if shiftFlag == 1
    [hei, wid, ~] = size( macSam );
    rrBioSam = zeros( hei, wid, 3, 'uint8' );
    [rh, rw, ~] = size( rBioSam );
    rrBioSam((hei-rh+1):end, (wid-rw+1):end, :) = rBioSam;
else
    rrBioSam = rBioSam;
end
if verbose == 1
    figure;imshowpair( macSam, rrBioSam );
    figure;imshowpair( macSam, rMSSam );
end
%     figure;imshowpair( rMSSam, rrBioSam );
figure; imshow( rrBioSam ); hold on; h = imshow( rMSSam );
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
alphaMat = alphaVal;
set( h, 'AlphaData', alphaMat );
fHandle = gcf;
end
