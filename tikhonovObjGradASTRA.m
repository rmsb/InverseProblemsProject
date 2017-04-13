% tikhonovObjGradASTRA.m

% We could propably do with less arguments, but lets start with simply
% takin everything
function result = tikhonovObjGradASTRA(f, m, fstar, alpha, ...
                                         proj_geom, vol_geom, proj_id, ...
                                         N, numAngles)
                                     
    % Sinogram row number
    nRow = length(m)/numAngles;

    % We want: 2*(A'*(A*f - m) + alpha*(f - fstar));
    % Replace A* with foward projection and A'* with
    % back projection.
    
    % Calculate FP = A*f
    FP_id = astra_mex_data2d('create', '-sino', proj_geom, 0);
    f_id = astra_mex_data2d('create', '-vol', vol_geom, reshape(f, N, N));
    cfgFP = astra_struct('FP');
    cfgFP.ProjectorId = proj_id;
    cfgFP.ProjectionDataId = FP_id;
    cfgFP.VolumeDataId = f_id;
    fp_id = astra_mex_algorithm('create', cfgFP);
    astra_mex_algorithm('run', fp_id); 
    FP = astra_mex_data2d('get', FP_id)';
    
    % Calculate P = A*f - m, in matrix form
    m = reshape(m, nRow, numAngles);
    P = FP - m;
    
    % Calculate BP = A'(A*f - m)
    P_id = astra_mex_data2d('create', '-sino', proj_geom, P');
    BP_id = astra_mex_data2d('create', '-vol', vol_geom, 0);
    cfgBP = astra_struct('BP');
    cfgBP.ProjectorId = proj_id;
    cfgBP.ProjectionDataId = P_id;
    cfgBP.ReconstructionDataId = BP_id;
    bp_id = astra_mex_algorithm('create', cfgBP);
    astra_mex_algorithm('run', bp_id);
    BP = astra_mex_data2d('get', BP_id);
    
    % Calculate the rest: 2*(BP + alpha*(f - fstar)
    result = 2*(BP(:) + alpha*(f - fstar));
    
end