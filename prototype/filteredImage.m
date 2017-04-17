function filtered = filteredImage(recon, radius, cutoff) 
% FILTEREDIMAGE applies predefined filters to a given pixel image. 

% Project values to the non-negative half. 
filtered = max(recon, 0);

% Apply circular median filter.

% This just creates a matrix with ones in a circular shape
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
