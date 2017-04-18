%% Reconstructions for 20-angles with 128x128 resolution.
% This file does the following.
% - Generates the measurement matrix.
% - Sparsifies the measurement sinogram.
% - Constructs the median-filtered target reconstruction. 
% - Computes the reconstructions with various alpha.
% - Produces the L-curve. 

%% Generate system matrix. 
sysmat = SystemMatrix(128, 20, 16);
save('data/SystemMatrix','sysmat');

%% Sparsify the sinogram. 
[sparseSino, sparseRecon] = createSparseSinogram(sinogram, sysmat);
 
%% Construct median-filtered target reconstruction. 
% TO-DO

%% Compute reconstructions with various values of alpha.
% TO-DO
fstar = zeros(sysmat.N);
alphas = ones(1);
[dataPenalty, regularizationPenalty] = computeReconstructions(...
    sysmat, sparseSino, fstar, alphas);

%% Produce the L-curve. 
figure();
plot(dataPenalty, regularizationPenalty);

