% Tikhonov_Lcurve_quadprog.m

%% Set the parameters

% Set the side length of the image, the downscale factor
% (from 2048) and the number of angles
downscaleFactor = 64;
N = 2048/downscaleFactor;
numAngles = 20;

% Use fstar = zero(N^2, 1) (Classical Tikhonov)? 1 for yes
fstarIsZero = 1;

% If fstarIsZero == 0, set the radius for the median filter 
% and the cutoff treshold
radius = 1;
cutoff = 0.15;

% Set the regularization parameter. For plotting nAlpha is set as a
% product.
xAlpha = 5;
yAlpha = 5;
nAlpha = xAlpha * yAlpha;
alpha = logspace(-2, 3, nAlpha);

% Options for quadprog
maxIterations = 300; % Default 200
optimalityTolerance = 10^-14; % Default 10^-8
constraintTolerance = 10^-12; % Default 10^-8
stepTolerance = 10^-14; % Default 10^-12


%% Load the data
loadcommand = ['load data\Sino_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' sinoShifted'];
eval(loadcommand);

loadcommand = ['load data\A_downsample', num2str(downscaleFactor), '_N', ...
               num2str(N), '_Ang', num2str(numAngles), ' A'];
eval(loadcommand);

loadcommand = ['load data\Recon_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' reconLoResSparse'];
eval(loadcommand);


%% Set up m, fstar and the functions for the PBB

% Drop m to vector
m = sinoShifted(:);

% Set the priori reconstruction
if fstarIsZero == 1
    fstar = zeros(N^2, 1);
else
    fstar = ellipticCut2(reconLoResSparse, [0.5, 0.51], 0.37, 0.27);
    fstar = applyFilters(fstar, radius, cutoff);
    
    figure(3)
    clf
    imagesc(fstar)
    axis image
    colorbar
    drawnow
  
    fstar = fstar(:);
end

% To speed up things alculate A'*A and A'*m. We seem to have enough
% memory for this.
AtA = A'*A;
Atm = A'*m;

% Define lower bound for f
lb = zeros(N^2, 1);
          

%% Quadprog
             
% Initialize the values for the L-curve
LCurveX = NaN(1, nAlpha);
LCurveY = NaN(1, nAlpha);

% Clear figure 1
figure(1)
clf
drawnow

% Compute the reconstruction for the alpha
for iii = 1:nAlpha
    % Follow what is happening
    disp([iii nAlpha])
    disp(alpha(iii))
    
    % Create H and h for specific alpha
    H = 2*(AtA + alpha(iii)*speye(N^2));
    h = -2*(Atm + alpha(iii)*fstar); 
    
    % Quadprog
    options = optimoptions(@quadprog, ...
                           'MaxIterations', maxIterations, ...
                           'OptimalityTolerance', optimalityTolerance, ...
                           'ConstraintTolerance', constraintTolerance, ...
                           'StepTolerance', stepTolerance);
    f = quadprog(H, h, [], [], [], [], lb, [], [], options);
    
    % Calculate the L-curve values
    LCurveX(iii) = log(norm(A*f - m));
    LCurveY(iii) = log(norm(f - fstar));
    
    % Plot the reconstructions
    figure(1)
    subplot(yAlpha, xAlpha, iii)
    imagesc(reshape(f, N, N))
    title(num2str(alpha(iii)))
    axis image
    axis off
    colormap default
    colorbar
    
    % Plot the Lcurve
    figure(2)
    clf
    plot(LCurveX, LCurveY)
    hold on
    plot(LCurveX, LCurveY, 'r.', 'MarkerSize', 14)
    axis equal
    drawnow
end
