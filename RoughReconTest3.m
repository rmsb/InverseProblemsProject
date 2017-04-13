% RoughReconTest3.m

%% Set the parameters

% Side length and number of angles
downscaleFactor = 64;
N = 2048/downscaleFactor;
numAngles = 20;

% Radius of meadian filter and lowpass cutoff treshold
nRad = 5;
radius = linspace(0.5, 2, nRad); % 10-25 for 512
nCut = 6;
lowPassCutoff = linspace(0.1, 0.2, nCut); % 10-25 for 512

%% Load the recon
loadcommand = ['load Recon_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' reconLoResSparse'];
eval(loadcommand);

%% Apply filters

% Set stuff outside a suitable ellips zero
recon = ellipticCut2(reconLoResSparse, [0.5, 0.51], 0.37, 0.27);

% Apply the rest of the filters
figure(1)
clf
drawnow
for iii = 1:nRad
    for jjj = 1:nCut        
        % Plot
        figure(1)
        subplot(nRad, nCut, (iii - 1)*nCut + jjj)
        imagesc(applyFilters(recon, radius(iii), lowPassCutoff(jjj)))
        axis image
        axis off
        title(['rad: ', num2str(radius(iii)), ', cut: ', num2str(lowPassCutoff(jjj))])
        drawnow
    end
end

