%------------------------------------------------------------------------
% This file is part of the
% All Scale Tomographic Reconstruction Antwerp Toolbox ("ASTRA-Toolbox")
%
% Copyright: iMinds-Vision Lab, University of Antwerp
% License: Open Source under GPLv3
% Contact: mailto:astra@ua.ac.be
% Website: http://astra.ua.ac.be
%------------------------------------------------------------------------
N = 32; % Edge length of reconstruction grid
f = 64; % Scaling factor for sinogram detector pixels
vol_geom = astra_create_vol_geom(N,N);

pixelsize = 0.050;
Dsd = 547.8; % source to detector
b = 204.3; 
x = 250; 

% source to origin
Dss = b + x; 
% Origin to detector
Otd = Dsd - Dss;
% Geometric magnification factor
M = Dsd / Dss; 
% Effective pixel size
effpixel = pixelsize / M; 
% Distance between source and center of rotation in pixels
D = Dss / effpixel; 
% Define the angles
% angInt = 10;
% angles = (1:angInt:360)*pi/180;
numAngles = 20;
angles = (1:360/numAngles:360)*pi/180;
% Construct the measurement geometry
proj_geom = astra_create_proj_geom('fanflat', M, 2240/f, angles, Dss/(effpixel*f), Otd/(effpixel*f));

proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);
% Generate the projection matrix for this projection model.
% This creates a matrix W where entry w_{i,j} corresponds to the
% contribution of volume element j to detector element i.
matrix_id = astra_mex_projector('matrix', proj_id);

% Get the projection matrix as a Matlab sparse matrix.
A = astra_mex_matrix('get', matrix_id);

% Save the system matrix A as a file
eval(['save A_downsample',num2str(f),'_N',num2str(N),'_Ang',num2str(numAngles), ' A'])

% % Test the reconstruction using lsqr
% P = phantom(N);
% sin = A*P(:);
% y = lsqr(A, sin(:), 1e-5, 100);
% rec = reshape(y, N,N);
% 
% figure;
% imagesc(reshape(sin,numel(sin)/numel(angles),numel(angles)));colormap gray;axis square
% figure;
% subplot(1,2,1)
% imagesc(P);colormap gray;axis square
% subplot(1,2,2)
% imagesc(rec);colormap gray;axis square

% Because Matlab's matrices are stored transposed in memory compared to C++,
% reshaping them to a vector doesn't give the right ordering for multiplication
% with W. We have to take the transpose of the input and output to get the same
% results (up to numerical noise) as using the toolbox directly.

% Each row of the projection matrix corresponds to a detector element.
% Detector t for angle p is for row 1 + t + p*proj_geom.DetectorCount.
% Each column corresponds to a volume pixel.
% Pixel (x,y) corresponds to column 1 + x + y*vol_geom.GridColCount.

astra_mex_projector('delete', proj_id);
astra_mex_matrix('delete', matrix_id);
