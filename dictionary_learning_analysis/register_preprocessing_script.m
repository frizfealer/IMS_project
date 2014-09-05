%% register_sImg_mImg_and_bImg_mImg
%register spectrum image to machine image
%and register bio image to machine image
%projectDir = the directory where the experiment is in.
%output: registeredInfo.mat
%% get the mass image with the most largest area
projectDir = '';
cd projectDir;
h = figure; imagesc( BlkDS.indMap );
axis off;
export_fig MSLoc.jpg
close;

[ sFileN, sFileP, ~] = uigetfile('*.jpg', 'Select the specturm location file' );
[ mFileN, mFileP, ~] = uigetfile('*.jpg', 'Select the machine image file' );
macSam = imread( strcat( mFileP, mFileN ) );
MSSam = imread( strcat( sFileP, sFileN ) );
fprintf( 'Arbitrary choosing points in machine image, because the correct information is in *.mis file\n' );
fprintf( 'save points as base_points and input_points\n' );
cpselect( MSSam, macSam );
fprintf( 'remember to change base_points by the location in *.mis file!!\n\n' );
%set the machine base_pnts to the points listed in *.mis file
mytForm = fitgeotrans( input_points, base_points, 'projective' );

%% testing the result...
mapMSSam = imwarp( MSSam, mytForm );
figure;imshowpair( macSam, mapMSSam );

%% testing the result with shifting...
fprintf( 'computing the shift of HEIGHT and WIDTH manually...\n' );
fprintf( 'HEIGHT is the y distance between the machine picture upper border and the plate''s upper border.\n' );
fprintf( 'WIDTH is the x distance between the machine picture left border and the plate''s left border.\n' );
shiftHEI = 11;
shiftWID = 44;
%loading the machine picture without margin
fprintf( 'edit the machine image file and save as fileName_2.jpg\n' );
[ mFileN2, mFileP2, ~] = uigetfile('*.jpg', 'Select the machine image file without margin' );
macSam = imread( strcat( mFileP2, mFileN2 ) );
[hei, wid, ~] = size( mapMSSam );
MMapMSSam = zeros( hei, wid, 3, 'uint8' );
MMapMSSam(1:(hei-shiftHEI), 1:(wid-shiftWID), :) = mapMSSam((shiftHEI+1):hei, (shiftWID+1):wid, :);
figure;imshowpair( macSam, MMapMSSam );

%% register bio image to machine image
[ sFileN, sFileP, ~] = uigetfile('*.jpg', 'Select the biological image file' );
[ mFileN, mFileP, ~] = uigetfile('*.jpg', 'Select the machine image file without margin' );
macSam = imread( strcat( mFileP, mFileN ) );
bioSam = imread( strcat( sFileP, sFileN ) );
cpselect( bioSam, macSam );
%save input_points1 and base_points1
mytForm2 = fitgeotrans( input_points1, base_points1, 'projective' );

%% testing the result...
mapbioSam = imwarp( bioSam, mytForm2 );
figure; imshow( macSam ); hold on; h = imshow( mapbioSam );alphaMat = 0.1;
set( h, 'AlphaData', alphaMat );

save('registeredInfo.mat', 'mytForm', 'mytForm2', 'shiftHEI', 'shiftWID' );