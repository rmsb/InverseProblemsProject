%% Load data matrices. 
filteredPriorRecon = load('data/Recon_128_prior.mat');
zeroPriorRecon = load('data/Recon_128_zeroprior.mat');
filteredPrior = load('data/InitialGuess');
fbpSparseRecon = load('prototype/data/SparseData');

%% Define color map axes. 

reconVals = [...
    filteredPriorRecon.recon(:); ...
    zeroPriorRecon.recon(:); ...
    filteredPrior.filteredRecon(:)];

colorRange = [min(reconVals), max(reconVals)];

%% Gather reconstructions and names into structs. 

recons = {filteredPriorRecon.recon;...
    zeroPriorRecon.recon;...
    filteredPrior.filteredRecon;...
    fbpSparseRecon.sparseRecon};

filenames = {'ReconFilteredPrior',...
    'ReconZeroPrior',...
    'FilteredPrior', ...
    'ReconSparseFBP'};

%% Plot reconstructions with titles to separate figures.

for i = 1:size(recons,1)
  currecon = recons{i};
  figure();
  imagesc(currecon);
  axis square;
  caxis(colorRange);
  xticks([]);
  yticks([]);
  filepath = strjoin({'images/', filenames{i}},'');
  saveas(gcf, filepath, 'epsc');
end
