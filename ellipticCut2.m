% ellipticCut2.m

function filtered = ellipticCut2(recon, centerPoint, majorAxis, minorAxis) 

% Image side length
N = size(recon, 1);

% Blah blah
[x, y] = meshgrid(linspace(0, 1, N), linspace(0, 1, N));

% Set points outside ellipse zero
filtered = recon;
filtered((x - centerPoint(1)).^2 / majorAxis^2 + ...
         (y - centerPoint(2)).^2 / minorAxis^2 > 1) = 0;

end