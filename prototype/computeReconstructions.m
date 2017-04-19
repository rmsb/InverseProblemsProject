function [dataPenalty, regularizationPenalty] = computeReconstructions(...
    A, sinogram, fstar, alphas)
% COMPUTECONSTRUCTIONS Regularized Tikhonov reconstructions for alphas.
%
% Computes argmin{|Af - m| + alpha * |f - fstar|} for each alpha. 
% For each reconstruction f, the log-penalty terms are calculated and
% returned in vectors. 
% 
% A - SystemMatrix object of tomography setting.
% sinogram - Measurement sinogram that matches the dimensions of A. 
% fstar - The regularization  model. 
% 

if isa(A, 'SystemMatrix')
    A= A.Matrix;
else
    error('Input argument A has to be an object of class SystemMatrix.');
end

m = sinogram(:);
MAXITER = 1000;
algorithm = ReconstructionAlgorithm(fstar, A, m, MAXITER);
nAlpha = length(alphas);
dataPenalty = NaN(1, nAlpha);
regularizationPenalty = NaN(1, nAlpha);

% Compute the reconstruction for the alpha
iteration = 0;
for iii = 1:nAlpha
    fprintf('Iteration: %d.\n', iteration);
    recon = algorithm.computeReconstruction(alphas(iii));
    dataPenalty(iii) = log(norm(A*recon - m));
    regularizationPenalty(iii) = log(norm(recon - fstar(:)));
    iteration = iteration + 1;
end

end

