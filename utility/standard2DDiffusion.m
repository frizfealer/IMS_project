function [ outW ] = standard2DDiffusion( inW, itNum, dt )
%UNTITLED3 Summary of this function goes here
% Solves the diffusion equation I_t = I_xx + I_yy
outW = inW;

% initialize the values (i+1,j), (i-1,j), (i,j+1), (i,j-1)
Ixp = zeros( size( outW ) );
Ixm = zeros( size( outW ) );
Iyp = zeros( size( outW ) );
Iym = zeros( size( outW ) );
% Ixpyp = zeros( size( outW ) );
% Ixmym = zeros( size( outW ) );
% Ixpym = zeros( size( outW ) );
% Ixmyp = zeros( size( outW ) );

iter = itNum; % number of iterations
% dt = 0.25;   % timestep

useNeumann = false; % use Neumann boundary conditions, otherwise Dirichlet

% iterate
for iI=0:iter

%   % display the image every 100 iteration steps
  if ( mod(iI,10)==0 )
%     fprintf( '%d\n', iI );
%     clf, imagesc( outW )
%     colormap( gray ), axis image, axis off
%     waitforbuttonpress
  end
  
  % compute the neighboring values 
  Ixp(1:end-1,:) = outW(2:end,:); Ixm(2:end,:) = outW(1:end-1,:);
  Iyp(:,1:end-1) = outW(:,2:end); Iym(:,2:end) = outW(:,1:end-1);
  
%   % for geometric equation neighboring values
%   Ixpyp(1:end-1, 1:end-1) = outW(2:end, 2:end);
%   Ixmym(2:end, 2:end) = outW(1:end-1, 1:end-1);
%   Ixpym(1:end-1,2:end) = outW(2:end, 1:end-1);
%   Ixmyp(2:end, 1:end-1) = outW(1:end-1, 2:end );
  
  % use Neumann boundary conditions (zero derivative at boundary)
  Ixp(end,:) = outW(end,:); Ixm(1,:) = outW(1,:);
  Iyp(:,end) = outW(:,end); Iym(:,1) = outW(:,1);
%   
%   % for geometric equation neighboring values
%   Ixpyp(end,:)=outW(end,:); Ixpyp(:,end)=outW(:,end);
%   Ixmym(1,:)=outW(1,:); Ixmym(:,1)=outW(:,1);
%   Ixpym(end,:)=outW(end,:); Ixpym(:,1)=outW(:,1);
%   Ixmyp(1,:)=outW(1,:); Ixmyp(:,end)=outW(:,end);

  % compute the derivatives
  Ixx = Ixp-2*outW+Ixm;
  Iyy = Iyp-2*outW+Iym;
%   Ixy = (Ixpyp + Ixmym - Ixpym -Ixmyp)./ 4;
%   Ix1 = Ixp - outW;
%   Iy1 = Iyp - outW;
  %ItHistMatrix(:,:,iI + 1) = Ixx + Iyy;
  % Geometric form
%   Ig = ( ( Iy1 .^2 ) .* Ixx - 2*Ixy .* Ix1 .* Iy1 + ( Ix1 .^2 ) .* Iyy ) ./ ( ( Ix1 .^2 ) + ( Iy1 .^2 ) + 1e-3 );
  % compute the new value by averaging
  Ikp1 = outW + dt*( Ixx + Iyy );
%   Ikp1 = outW + dt*( Ig );
  if ( useNeumann )
    outW = Ikp1;
  else  % using Dirichlet boundary conditions (ie., fix the boundary)
    outW(2:end-1,2:end-1) = Ikp1(2:end-1,2:end-1);
  end
  
end
outW = outW(:);
end

