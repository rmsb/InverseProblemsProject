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
load data/SparseData sparseSino
load data/SystemMatrix sysmat
nAlpha = 50;
base = 1.2;
aMin = 10^-4;
aMax = 500;
alphas = base.^linspace(log(aMin)/log(base), log(aMax)/log(base), nAlpha);
fstar = zeros(sysmat.N);
[dataPenalty, regularizationPenalty] = computeReconstructions(...
    sysmat, sparseSino, fstar, alphas);
save('data/Penalties', 'dataPenalty', 'regularizationPenalty', ...
    'alphas', 'fstar');

%% Produce the L-curve and select value of alpha. 
penalties = load('data/Penalties');
alphaIndex = 35;

figure();
plot(penalties.dataPenalty, penalties.regularizationPenalty, '-');
hold on;
plot(penalties.dataPenalty(alphaIndex), ...
    penalties.regularizationPenalty(alphaIndex), 'or');
alpha = penalties.alphas(alphaIndex);
hold off;

%% Create the optimal reconstruction.
load data/SystemMatrix sysmat
load data/SparseData sparseSino
penalties = load('data/Penalties');

MAXITER = 1000;
algorithm = ReconstructionAlgorithm(penalties.fstar, sysmat.Matrix, ...
    sparseSino, MAXITER);
optimalRecon = algorithm.computeReconstruction(alpha);
optimalRecon = reshape(optimalRecon, sysmat.N*[1,1]);
save('data/OptimalRecon', 'optimalRecon', 'alpha');

%% Plot the optimal reconstruction. 
optimRecon = load('data/OptimalRecon');
imagesc(optimRecon.optimalRecon);
colorbar;
