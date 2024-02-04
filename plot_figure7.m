function plot_figure7(true, error_l_pca, error_l_gKDR, error_c,stdg_l_pca,stdg_l_gKDR, stdg_c,  predicting, units)
% Plots true values and root mean square errors (RMSE) for predictions.
%
% Parameters:
% true:       matrix with true values
% error_l:    error values for linked emulator
% error_c:    error values for composite emulator
% predicting: the type of value being predicted (e.g., 'Pressure')
% units:      the units for the values being plotted (e.g., '(psi)')

% Create a new figure with default font size 18 for axes
figure('DefaultAxesFontSize', 16)
set(gcf, 'Position', get(0, 'Screensize'), 'Visible', 'on');
%figure('DefaultAxesFontSize', 18)
tiledlayout(1,3);
nexttile
% First subplot: plotting true values
plot(0:length(true) - 1, true');

ylabel(strcat(predicting, ' ', units), 'FontSize', 20) % y-axis label
xlabel('Depth (in)', 'FontSize', 20)                   % x-axis label
xlim([0, 100]);
% Set x-axis ticks and labels
xticks([0 10 20 30 40 50 60 70 80 90 100])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'})

view([90 90])                                          % Rotate plot for a horizontal view
pbaspect([2 1 1]) % Set plot aspect ratio

% Maximize figure on screen
%figure('DefaultAxesFontSize', 18)
%set(gcf, 'Position', get(0, 'Screensize'), 'Visible', 'on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
% Second subplot: plotting RMSE for the two emulators
% Calculate and plot RMSE for composite emulator
plot(0:length(sqrt(mean(error_c.^2))) - 1, sqrt(mean(error_c.^2)), 'LineWidth', 3);
hold on
% Calculate and plot RMSE for linked emulator
plot(0:length(sqrt(mean(error_l_pca.^2))) - 1, sqrt(mean(error_l_pca.^2)), 'Color', 'k', 'LineWidth', 3,'LineStyle','--');
plot(0:length(sqrt(mean(error_l_gKDR.^2))) - 1, sqrt(mean(error_l_gKDR.^2)), 'Color', 'r', 'LineWidth', 3,'LineStyle',':');

ylabel(strcat('RMSE (psi)'), 'FontSize', 20)   % y-axis label
xlabel('Depth (in)', 'FontSize', 20)                  % x-axis label
%legend('PPCE', 'PPLE', 'Location', 'best', 'FontSize', 20) % Add legend
ylim([0 6e6])
xlim([0, 100])
% Set x-axis ticks and labels
xticks([0 10 20 30 40 50 60 70 80 90 100])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'})
pbaspect([2 1 1]) % Set plot aspect ratio
view([90 90]) % Rotate plot for a horizontal view

% %plot credible intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile
av_LCI_c = 2*1.96*mean(stdg_c);
av_LCI_l_pca =  2*1.96*mean(stdg_l_pca);
av_LCI_l_gKDR =  2*1.96*mean(stdg_l_gKDR);
% Plotting for composite emulator
plot(0:length(av_LCI_c) - 1,av_LCI_c, 'LineWidth', 3);
hold on
% Plotting for linked emulator
plot(0:length(av_LCI_l_pca) - 1, av_LCI_l_pca, 'Color', 'k', 'LineWidth', 3,'LineStyle','--');
plot(0:length(av_LCI_l_gKDR) - 1, av_LCI_l_gKDR, 'Color', 'r', 'LineWidth', 3,'LineStyle',':');

ylabel(strcat('L_{CI} (psi)'), 'FontSize', 20)   % y-axis label
xlabel('Depth (in)', 'FontSize', 20)                  % x-axis label
%legend('PPCE', 'PPLE', 'Location', 'best', 'FontSize', 20) % Add legend

xlim([0, 100])
% Set x-axis ticks and labels
xticks([0 10 20 30 40 50 60 70 80 90 100])
pbaspect([2 1 1])
xticklabels({'0', '10', '20', '30', '40', '50', '60', '70', '80', '90', '100'})

view([90 90]) % Rotate plot for a horizontal view


end