% CreateData.m

% - Downscale the sinogram, create a sparse version and save.
% - Create ifanbeam reconstruction and save.
% - Shift the columns of the sinogram, so the astra matrix reconstruction
%   has the right orientation. (We could also shift the rows
%   of the astra matrix.) This has to be done after the recons
%   are donw with ifanbeam.

%% Set the parameters

% Same parameters as in BuildMatrixASTRA.m.
% I think N = 2048/downscaleFactor works quite well.
downscaleFactor = 8; % f in BuildMatrixASTRA.m
N = 2048/downscaleFactor;
numAngles = 20;

%% Load data
load LotusSinogram

% Load measurement matrix A for comparing shifted sinogram and A*recon
loadcommand = ['load A_downsample', num2str(downscaleFactor), '_N', ...
               num2str(N), '_Ang', num2str(numAngles), ' A'];
eval(loadcommand);


%% Downscale the sinogram and create a sparse version

figure(1);
clf
subplot(3, 1, 1)
imagesc(sinogram);
colorbar

% Sum downscaleFactor number of rows
[rows, cols] = size(sinogram);
scalingMatrix = kron(speye(rows/downscaleFactor*cols), ones(1, downscaleFactor));
sinoLoRes = scalingMatrix * sinogram(:);
sinoLoRes = reshape(sinoLoRes, rows/downscaleFactor, cols);

figure(1);
subplot(3, 1, 2)
imagesc(sinoLoRes);
colorbar

% Create the sparse version
sinoLoResSparse = sinoLoRes(:, 1:360/numAngles:360);

figure(1)
subplot(3, 1, 3)
imagesc(sinoLoResSparse);
colorbar


%% Specify measurement geometry

% Size of detector pixel in millimeters
pixelsize = 0.050;

% Distance from X-ray source to X-ray detector in mm
% Always has the same constant value
Dsd = 547.8;

% Distance from X-ray source to home point of translation stage in mm
% Always has the same constant value
b = 204.3;

% Distance from home point of translation stage to position of stage in mm
% This is a variable
x = 250;

% Distance from X-ray source to center of rotation in mm
Dss = b + x;

% Geometric magnification factor
M = Dsd / Dss;

% Effective pixel size that takes into account the geometric magnification
effpixel = pixelsize / M;

% Distance between X-ray source and center of rotation expressed in pixels
D   = Dss / effpixel;


%% Create the ifanbeam reconstructions

recon = ifanbeam(sinogram, D, ...
    'FanRotationIncrement', 1, ...
    'FanSensorGeometry', 'line', ...
    'OutputSize', 2048);

reconLoRes = ifanbeam(sinoLoRes, D/downscaleFactor, ...
    'FanRotationIncrement', 1, ...
    'FanSensorGeometry', 'line', ...
    'OutputSize', N);

reconLoResSparse = ifanbeam(sinoLoResSparse, D/downscaleFactor, ...
    'FanRotationIncrement', (360/numAngles), ...
    'FanSensorGeometry', 'line', 'OutputSize', N);

figure(2)
clf
subplot(1, 3, 1)
imagesc(recon)
axis image
axis off
colorbar


subplot(1, 3, 2)
imagesc(reconLoRes)
axis image
axis off
colorbar

subplot(1, 3, 3)
imagesc(reconLoResSparse)
axis image
axis off
colorbar


%% Shift the rows so..
% In the usual radon 0 angle rays go up and we move counter clockwise.
% It seems here zero angle rays go left and we move clockwise.
sinoShifted = sinoLoResSparse;

% Rotate 90 degrees
% Index of the last angle in the first quarter. Seems to work with
% nice angle numbers (such as 20) and I don't want to think about
% correct rounding in ugly cases..
quarterInd = numAngles/4 + 1;
sinoShifted = [sinoShifted(:, quarterInd + 1:end), ...
                   sinoShifted(:, 1:quarterInd)];

% Change rotation direction
sinoShifted = fliplr(sinoShifted);

% Create sinogram with measurement matrix and FBP recon.
% This rough shape is what we are trying to achieve by shifting.
sinoARecon = reshape(A*reconLoRes(:), size(A, 1)/numAngles, numAngles);

% Compare visually
figure(3)
subplot(2, 1, 1)
imagesc(sinoARecon)
colorbar

subplot(2, 1, 2)
imagesc(sinoShifted)
colorbar


%% Save in the same format as the measurement matrix
% Save the sparse sinogram
savecommand = ['save Sino_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' sinoShifted'];
eval(savecommand);

% Save the sparse reconstruction
savecommand = ['save Recon_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' reconLoResSparse'];
eval(savecommand);
