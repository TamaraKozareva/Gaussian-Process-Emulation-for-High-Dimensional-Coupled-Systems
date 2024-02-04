function [error_r, true_r, predg_r, stdg_r] = save_data(predg, true_pres, stdg)
% Calculates and returns prediction errors and retains relevant data.
%
% Parameters:
% predg:      predicted values.
% true_pres:  true values.
% stdg:       standard deviations of the predicted values.
%
% Returns:
% error_r:    errors between predicted and true values.
% true_r:     true values.
% predg_r:    predicted values.
% stdg_r:     standard deviations of the predicted values.

% Calculate the error between predicted and true values
error_r = predg - true_pres;

% Store the other provided values directly into the return variables
true_r = true_pres;
predg_r = predg;
stdg_r = stdg;

end