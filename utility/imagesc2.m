function h = imagesc2 ( img_data )
% a wrapper for imagesc, with some formatting going on for nans

% plotting data. Removing and scaling axes (this is for image plotting)
h = imagesc(img_data);
% axis image off

% setting alpha values
if ndims( img_data ) == 2
  set(h, 'AlphaData', img_data~=0)
elseif ndims( img_data ) == 3
  set(h, 'AlphaData', ~isnan(img_data(:, :, 1)))
end
  set(gca, 'Color', [0, 0, 0] );
if nargout < 1
  clear h
end