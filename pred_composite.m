function [mean,stdg]= pred_composite(xxall, newu, gall, model)
% Predicts the mean and standard deviation for the given data and model using composite emulator.
%
% Parameters:
% xxall:  Training input data.
% newu:   Test input data for prediction.
% gall:   Output data corresponding to xxall.
% model:  Structure containing parameters and settings of the model.
%
% Returns:
% mean:   Predicted mean.
% stdg:   Predicted standard deviation.

% Dimension information
dimz = size(gall,2); % Dimension of the output
dimy = 0;            % Initializing dimension of another possible output (unused)

% Calculate the Ru matrix using the kernel function
Ru = kernel_PE(xxall, xxall, model.alpha(1), model.range_par);
% Adding nugget term for regularization/stability
Ru = Ru + eye(size(Ru)) * (model.nugget);

% Invert the Ru matrix
invRu = inv(Ru);

% Calculate the r vector
r = kernel_PE(xxall, newu, model.alpha(1), model.range_par);

% Construct the Hxdg matrix for the design points
xdg = xxall;
Hxdg = [ones(size(xxall,1),1) xdg];

% Construct the Hxdg matrix for prediction points
Hxdg_star = [ones(size(newu,1),1) newu];

% Calculate regression coefficients thetas
thetas = (Hxdg' / Ru * Hxdg) \ (Hxdg' / Ru * gall);

% Predict the mean
mean = Hxdg_star * thetas + r' / Ru * (gall - Hxdg * thetas);

% Calculate sigma squared
sigma_sqred = diag(sigma_sqr(Hxdg, invRu, gall, thetas, size(gall,1), size(thetas,1)));

% Compute the variance for predictions
var = sigma_sqred' .* (1 + model.nugget - r' / Ru * r + (Hxdg_star' - Hxdg' / Ru * r)' / (Hxdg' / Ru * Hxdg) * (Hxdg_star' - Hxdg' / Ru * r));

% Calculate standard deviation from variance
stdg = sqrt(var);

end