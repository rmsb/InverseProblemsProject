%% Plots the LCurves for the two different reconstructions.

filterCurve = load('data/LCurve_128_prior');
zeroCurve = load('data/LCurve_128_zeroprior');
filterCorner = load('data/Recon_128_prior'); 
zeroCorner = load('data/Recon_128_zeroprior');

%% Prepare structs. 

LCurveX = {filterCurve.LCurveX; zeroCurve.LCurveX};
LCurveY = {filterCurve.LCurveY; zeroCurve.LCurveY};
CornerX = {filterCorner.LCurveX; zeroCorner.LCurveX};
CornerY = {filterCorner.LCurveY; zeroCorner.LCurveY};
alpha = {filterCorner.alpha; zeroCorner.alpha};
filenames = {'LCurveFilteredPrior', 'LCurveZeroPrior'};

%% Make plots. 

for i = 1:size(LCurveX, 1)
    figure();
    plot(LCurveX{i}, LCurveY{i}, 'LineWidth', 5);
    hold on; 
    plot(CornerX{i}, CornerY{i}, 'r.', 'MarkerSize', 40);
    text(CornerX{i}, CornerY{i}, ['$\quad\:\alpha=', num2str(alpha{i}, 2), '$'], 'Interpreter', 'latex') 
    axis square;
    %xticks([]);
    %yticks([]);
    xlabel('$\log\left\Vert A\mathbf{f}-\mathbf{m}\right\Vert _2$', 'Interpreter', 'latex');
    ylabel('$\log\left\Vert \mathbf{f}-\mathbf{f}^\star\right\Vert _2$', 'Interpreter', 'latex');
    filepath = strjoin({'images/', filenames{i}}, '');
    saveas(gcf, filepath, 'epsc');
end
