% Tikhonov_Lcurve_PBB_ASTRA.m

%% Set the parameters

% Set the side length of the image, the downscale factor
% (from 2048) and the number of angles
downscaleFactor = 8;
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
xAlpha = 1;
yAlpha = 1;
nAlpha = xAlpha * yAlpha;
alpha = logspace(2, 2, nAlpha);

% Set maximum number of iterations in the PBB
MAXITER = 20;      


%% Load the data
loadcommand = ['load data\Sino_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' sinoShifted'];
eval(loadcommand);

loadcommand = ['load data\Recon_downsample', num2str(downscaleFactor), '_N', ...
                num2str(N), '_Ang', num2str(numAngles), ' reconLoResSparse'];
%eval(loadcommand);


%% Specify geometries, projector and data

% Create volume geometry
vol_geom = astra_create_vol_geom(N,N);

pixelsize = 0.050;
Dsd = 547.8; % source to detector
b = 204.3; % distance to "rails"
x = 250; % distance along th rails

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
angles = (1:360/numAngles:360)*pi/180;

% Construct the measurement geometry
proj_geom = astra_create_proj_geom('fanflat', M, 2240/downscaleFactor, ...
                                   angles, Dss/(effpixel*downscaleFactor), ...
                                   Otd/(effpixel*downscaleFactor));

% Create projector
proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

% Set the sinogram
sino_id = astra_mex_data2d('create', '-sino', proj_geom, sinoShifted');



%% Set up m, fstar and the functions for the PBB

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

% Drop sino as a vector
m = sinoShifted(:);
   
          
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
    
    gold   = tikhonovObjGradASTRA(fold, m, fstar, alpha(iii), ...
                                      proj_geom, vol_geom, proj_id, ...
                                      N, numAngles);
    obj(1) = tikhonovObjectiveASTRA(fold, m, fstar, alpha(iii), ...
                                    proj_geom, vol_geom, proj_id, ...
                                    N, numAngles);
    
    % Make the first iteration step. Theoretically, this step should satisfy
    % the Wolfe condition, see [J.Nocedal, Acta Numerica 1992].
    % We use simply a constant choice since it usually works well.
    % If there is a problem with convergence, try making t smaller.
    t = .000001;
    
    % Compute new iterate point fnew and gradient gnew at fnew
    fnew = max(fold - t*gold,0); % THE NON-NEGATIVITY PROJECTION of the first step
    gnew   = tikhonovObjGradASTRA(fnew, m, fstar, alpha(iii), ...
                                  proj_geom, vol_geom, proj_id, ...
                                  N, numAngles);
    

    % Iteration counter
    its = 2;
    
    % Record value of objective function at the new point
    obj(its) = tikhonovObjectiveASTRA(fnew, m, fstar, alpha(iii), ...
                                 proj_geom, vol_geom, proj_id, ...
                                 N, numAngles);
    
    % Follow if objective function has minimal value so far
    % and save the f that produces the minimal value. The initial
    % values don't really matter.
    objMin = obj(its);
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
        
        gnew   = tikhonovObjGradASTRA(fnew, m, fstar, alpha(iii), ...
                                          proj_geom, vol_geom, proj_id, ...
                                          N, numAngles);
        obj(its) = tikhonovObjectiveASTRA(fnew, m, fstar, alpha(iii), ...
                                          proj_geom, vol_geom, proj_id, ...
                                          N, numAngles);
        
        % Follow what is happening
        if (obj(its) < objMin)
            fMin = fnew;
            objMin = obj(its);
        end
        format short e
        
        % Print stuff
        disp(['Iteration ', num2str(its,'%4d'), ', objective function value ', ...
               num2str(obj(its), '%.3e'), ', difference: ', ... 
               num2str(obj(its) - obj(its - 1), '%.3e'), ...
               ', is smallest so far: ', num2str(objMin == obj(its))])
        
        % Plot the progress of the objective function
        figure(4)
        plot(obj)
        drawnow    
    end 
    
    % Pick the minimal iteration as recon
    recon = fMin;
    
    % Calculate FP = A*recon
    Arecon_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
    recon_id = astra_mex_data2d('create', '-vol', vol_geom, reshape(recon, N, N));
    cfgFP = astra_struct('FP');
    cfgFP.ProjectorId = proj_id;
    cfgFP.ProjectionDataId = Arecon_id;
    cfgFP.VolumeDataId = recon_id;
    fp_id = astra_mex_algorithm('create', cfgFP);
    astra_mex_algorithm('run', fp_id); 
    Arecon = astra_mex_data2d('get', Arecon_id)';
    
    % Calculate the L-curve values
    LCurveX(iii) = log(norm(Arecon(:) - m));
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
