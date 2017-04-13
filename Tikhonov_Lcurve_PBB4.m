% Tikhonov_Lcurve_PBB3.m

%% Set the parameters

% Set the side length of the image, the downscale factor
% (from 2048) and the number of angles
downscaleFactor = 16;
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
alpha = logspace(-3, 0, nAlpha);

% Set maximum number of iterations in the PBB
MAXITER = 20;      


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

% Create the objective function (Tikohonov penalty) and its gradient.
% I'm not sure if this particular use of brackets makes things faster of
% slower..
objective = @(f, A, AtA, m, Atm, fstar, alpha) ...
            f'*(AtA*f - 2*Atm + alpha*(f - 2*fstar));
objGradient = @(f, A, AtA, m, Atm, fstar, alpha) ...
              2*(AtA*f - Atm + alpha*(f - fstar));
   
          
%% The PBB
             
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
    
    % Optimization routine
    obj    = NaN(MAXITER + 1, 1);     % We will monitor the value of the objective function
    fold   = zeros(N^2, 1);    % Initial guess
    gold   = objGradient(fold, A, AtA, m, Atm, fstar, alpha(iii));
    obj(1) = objective(fold, A, AtA, m, Atm, fstar, alpha(iii));
    
    % Make the first iteration step. Theoretically, this step should satisfy
    % the Wolfe condition, see [J.Nocedal, Acta Numerica 1992].
    % We use simply a constant choice since it usually works well.
    % If there is a problem with convergence, try making t smaller.
    t = .000001;
    
    % Compute new iterate point fnew and gradient gnew at fnew
    fnew = max(fold - t*gold,0); % THE NON-NEGATIVITY PROJECTION of the first step
    gnew =  objGradient(fnew, A, AtA, m, Atm, fstar, alpha(iii));

    % Iteration counter
    its = 1;
    
    % Record value of objective function at the new point
    OFf        = objective(fnew, A, AtA, m, Atm, fstar, alpha(iii));
    obj(its+1) = OFf;
    
    % Follow if objective function has minimal value so far
    % and save the f that produces the minimal value. The initial
    % values don't really matter.
    objMin = OFf;
    fMin = fnew;

    % Barzilai and Borwein iterative minimization routine
    figure(4)
    clf
    while (its  < MAXITER)
        % Increment iteration counter
        its = its + 1;
        
        % Compute steplength alpha
        fdiff   = fnew - fold;
        gdiff   = gnew - gold;
        steplen = (fdiff'*fdiff)/(gdiff'*fdiff);
        
        % Update points, gradients and objective function value
        fold = fnew;
        gold = gnew;
        fnew = max(fnew - steplen*gnew,0);  % THE NON-NEGATIVITY PROJECTION
        gnew = objGradient(fnew, A, AtA, m, Atm, fstar, alpha(iii));
        OFf  = objective(fnew, A, AtA, m, Atm, fstar, alpha(iii));
        
        % Follow what is happening
        obj(its+1) = OFf;
        if (OFf < objMin)
            fMin = fnew;
            objMin = OFf;
        end
        format short e
        
        % Print stuff
        disp(['Iteration ', num2str(its,'%4d'), ', objective function value ', ...
               num2str(obj(its), '%.3e'), ', difference: ', ... 
               num2str(obj(its) - obj(its - 1), '%.3e'), ...
               ', is smallest so far: ', num2str(objMin == OFf)])
        
        % Plot the progress of the objective function
        figure(4)
        plot(obj)
        drawnow    
    end 
    
    % Pick the minimal iteration as recon
    recon = fMin;
    
    % Calculate the L-curve values
    LCurveX(iii) = log(norm(A*recon - m));
    LCurveY(iii) = log(norm(recon - fstar));
    
    % Plot the reconstructions
    figure(1)
    subplot(yAlpha, xAlpha, iii)
    imagesc(reshape(recon, N, N))
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
