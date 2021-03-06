% applyFilters.m

function filtered = applyFilters(recon, radius, cutoff) 

% Project to non-negative
filtered = max(recon, 0);

% Apply circular median filter
% This just creates a matrix with ones in a cricular shape
domain = fspecial('disk', radius) > 0;
% How many ones are in the matrix
nElements = numel(find(domain));
% The ordinal of the middle element
middle = round(0.5*nElements);
% Apply the median filter
filtered = ordfilt2(filtered, middle, domain);

% Apply the cutoff
filtered(filtered < max(filtered(:))*cutoff) = 0;

end