%% Reconstructions for 20-angles with 128x128 resolution.
% This file does the following.
% - Generates the measurement matrix.
% - Sparsifies the measurement sinogram.
% - Constructs the median-filtered target reconstruction. 
% - Computes the reconstructions with various alpha.
% - Produces the L-curve. 

%% Load sinogram. 
load data/LotusSinogram sinogram

%% Generate system matrix. 
sysmat = SystemMatrix(128, 20, 16);
save('data/SystemMatrix','sysmat');

%% Sparsify the sinogram. 
load data/SystemMatrix sysmat
[sparseSino, sparseRecon] = createSparseSinogram(sinogram, sysmat);
save('data/SparseData', 'sparseSino', 'sparseRecon'); 

%% Construct filtered target reconstruction. 
load data/SparseData sparseRecon
filteredRecon = ellipticCut(sparseRecon, [0.5, 0.51], 0.37, 0.27);
filteredRecon = filteredImage(filteredRecon, 3, 0.1);
save('data/InitialGuess', 'filteredRecon');

%% View filtered target reconstruction. 
load data/InitialGuess filteredRecon
imagesc(filteredRecon);

%% Compute reconstructions with various values of alpha.
% TO-DO
fstar = zeros(sysmat.N);
alphas = ones(1);
[dataPenalty, regularizationPenalty] = computeReconstructions(...
    sysmat, sparseSino, fstar, alphas);

%% Produce the L-curve. 
figure();
plot(dataPenalty, regularizationPenalty);

