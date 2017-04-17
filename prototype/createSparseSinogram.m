function [sparseSino, sparseRecon] = createSparseSinogram(sinogram, A)
% CREATESPARSESINOGRAM Sparsifies the given sinogram data to A's dimension. 
% 
% The function does the following.
% - Downscale the sinogram, create a sparse version and save.
% - Create ifanbeam reconstruction and save.
% - Shift the columns of the sinogram, so the astra matrix reconstruction
%   has the right orientation. (We could also shift the rows
%   of the astra matrix.) This has to be done after the recons
%   are done with ifanbeam.

if isa(A, 'SystemMatrix')
    N = A.N;
    numAngles = A.Angles;
    downscaleFactor = A.DownscaleFactor;
else
    error('Input argument A has to be an object of class SystemMatrix.');
end

sinoLoRes = downscaleMatrix(sinogram, downscaleFactor);
sinoLoResSparse = sinoLoRes(:, 1:360/numAngles:360);
sparseSino = shiftSinogram(sinoLoResSparse, numAngles);
sparseRecon = fanbeamReconstruction(sinoLoResSparse, downscaleFactor);
save(createFilePath(N, downscaleFactor, numAngles), 'sparseSino')
end

function pathstring = createFilePath(N, f, numAngles)
% createFilePath creates the path where the data matrix will be saved. 
components = {'data\SparseSinogram',num2str(f),'_N',...
    num2str(N),'_Ang',num2str(numAngles)};
pathstring  = strjoin(components, '');
end

function output = shiftSinogram(sinogram, numAngles)
% In the usual radon 0 angle rays go up and we move counter clockwise.
% It seems here zero angle rays go left and we move clockwise.

% Rotate 90 degrees
% Index of the last angle in the first quarter. Seems to work with
% nice angle numbers (such as 20) and I don't want to think about
% correct rounding in ugly cases..
quarterInd = numAngles/4 + 1;
sinoShifted = [sinogram(:, quarterInd + 1:end), ...
                   sinogram(:, 1:quarterInd)];
% Change rotation direction
output = fliplr(sinoShifted);
end

function output = downscaleMatrix(sinogram, downscaleFactor)
% Sum downscaleFactor number of rows
[rows, cols] = size(sinogram);
scalingMatrix = kron(speye(rows/downscaleFactor*cols), ...
    ones(1, downscaleFactor));
sinoLoRes = scalingMatrix * sinogram(:);
output = reshape(sinoLoRes, rows/downscaleFactor, cols);
end

function recon = fanbeamReconstruction(sinogram, downscaleFactor)
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

recon = ifanbeam(sinogram, D/downscaleFactor, ...
    'FanRotationIncrement', (360/numAngles), ...
    'FanSensorGeometry', 'line', 'OutputSize', N);
end
