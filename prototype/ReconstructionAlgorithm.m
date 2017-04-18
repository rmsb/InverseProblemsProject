classdef ReconstructionAlgorithm
    % RECONSTRUCTIONALGORITHM Projected Barzilai-Borwein algorithm object. 
    %
    % An object of this class can be used to calculate reconstructions with
    % the computeReconstruction() method. The system matrix, measurement 
    % data and fstar are properties of the object that are given in the
    % constructor. The computeReconstruction() method uses these property
    % values in calculating the reconstruction. 
    
    properties
        fstar
        A
        AtA
        Atm
        m
        MAXITER
    end
    
    methods
        
        function ra = ReconstructionAlgorithm(fstar, A, m, MAXITER)
            if nargin > 0
             ra.fstar = fstar(:);
             ra.A = A;
             ra.AtA = A'*A;
             ra.Atm = A'*m;
             ra.m = m(:);
             ra.MAXITER = MAXITER;
            end
        end
        
        function recon = computeReconstruction(obj, alpha)
            objective = NaN(obj.MAXITER + 1, 1);
            fold = zeros(numel(obj.fstar), 1);   
            gold = obj.objGradient(fold, alpha);
            objective(1) = obj.objective(fold, alpha);
            % Make the first iteration step. Theoretically, this step 
            % should satisfy the Wolfe condition, see [J.Nocedal, Acta 
            % Numerica 1992]. We use simply a constant choice since it 
            % usually works well. If there is a problem with convergence, 
            % try making t smaller.
            t = .000001;
            % Compute new iterate point fnew and gradient gnew at fnew
            % THE NON-NEGATIVITY PROJECTION of the first step
            fnew = max(fold - t*gold, 0);
            gnew = obj.objGradient(fnew, alpha);
            % Iteration counter
            its = 1;
            % Record value of objective function at the new point
            OFf = obj.objective(fnew, alpha);
            objective(its+1) = OFf;
            % Follow if objective function has minimal value so far
            % and save the f that produces the minimal value. The initial
            % values don't really matter.
            objMin = OFf;
            fMin = fnew;
            % Barzilai and Borwein iterative minimization routine
            while (its  < obj.MAXITER)
                % Increment iteration counter
                its = its + 1;
                % Compute steplength alpha
                fdiff = fnew - fold;
                gdiff = gnew - gold;
                steplen = (fdiff'*fdiff)/(gdiff'*fdiff);
                % Update points, gradients and objective function value
                fold = fnew;
                gold = gnew;
                % THE NON-NEGATIVITY PROJECTION
                fnew = max(fnew - steplen*gnew,0);  
                gnew = obj.objGradient(fnew, alpha);
                OFf = obj.objective(fnew, alpha);
                % Follow what is happening
                objective(its+1) = OFf;
                if (OFf < objMin)
                    fMin = fnew;
                    objMin = OFf;
                end
            end 
            % Pick the minimal iteration as recon
            recon = fMin;
        end
        
        function output = objective(obj, f, alpha)
            output = f'*(obj.AtA*f - 2*obj.Atm + alpha*(f - 2*obj.fstar));
        end
        
        function output = objGradient(obj, f, alpha)
            output = 2*(obj.AtA*f - obj.Atm + alpha*(f - obj.fstar));
        end
        
        function [] = printFstar(obj)
            disp(obj.fstar);
        end
    end
    
end

