classdef SystemMatrix
    %SYSTEMMATRIX Tomographic measurement matrix object. 
    %
    % Objects of this class create the system matrix for tomographical 
    % imaging using the fanbeam geometry of the device at the Physicum lab.
    % The constructor is given the properties of the setting and the system
    % matrix is constructed automatically using the ASTRA toolbox. 
    
    properties
        Matrix
        N
        Angles
        DownScaleFactor
    end
    
    methods
        function sm = SystemMatrix(N, Angles, DownScaleFactor)
            if nargin > 0
             sm.N = N;
             sm.Angles = Angles;
             sm.DownScaleFactor = DownScaleFactor;
             sm.Matrix = sm.calculateMatrix();
            end
        end
        
        function sysmat = calculateMatrix(obj)
            vol_geom = astra_create_vol_geom(obj.N, obj.N);
            pixelsize = 0.050;
            Dsd = 547.8; % source to detector
            b = 204.3; 
            x = 250; 
            % source to origin
            Dss = b + x; 
            % Origin to detector
            Otd = Dsd - Dss;
            % Geometric magnification factor
            M = Dsd / Dss; 
            % Effective pixel size
            effpixel = pixelsize / M; 
            % Define the angles
            angles = (1:360/obj.Angles:360)*pi/180;
            % Construct the measurement geometry
            f = obj.DownScaleFactor;
            proj_geom = astra_create_proj_geom('fanflat', M, 2240/f, ...
                angles, Dss/(effpixel*f), Otd/(effpixel*f));
            proj_id = astra_create_projector('line_fanflat', ...
                proj_geom, vol_geom);
            % Generate the projection matrix for this projection model.
            % This creates a matrix W where entry w_{i,j} corresponds to the
            % contribution of volume element j to detector element i.
            matrix_id = astra_mex_projector('matrix', proj_id);
            % Get the projection matrix as a Matlab sparse matrix.
            sysmat = astra_mex_matrix('get', matrix_id);
            % Because Matlab's matrices are stored transposed in memory compared to C++,
            % reshaping them to a vector doesn't give the right ordering for multiplication
            % with W. We have to take the transpose of the input and output to get the same
            % results (up to numerical noise) as using the toolbox directly.
            % Each row of the projection matrix corresponds to a detector element.
            % Detector t for angle p is for row 1 + t + p*proj_geom.DetectorCount.
            % Each column corresponds to a volume pixel.
            % Pixel (x,y) corresponds to column 1 + x + y*vol_geom.GridColCount.
            astra_mex_projector('delete', proj_id);
            astra_mex_matrix('delete', matrix_id);
        end
        
    end
    
end
