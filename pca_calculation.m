% pca_calculation function performs Principal Component Analysis on the given data 
% and returns the principal components based on the user-defined number of modes.
function [z, num_modes, total_var_z] = pca_calculation(data, user_modes)

    % Compute the mean of each column (feature) in the data.
    Xbar = mean(data, 1);
    
    % Center the data by subtracting the mean.
    X = (data - Xbar)';
    
    % Calculate the covariance matrix A.
    A = X * X';
    
    % Compute the Singular Value Decomposition (SVD) of the centered data.
    [E, D, V] = svd(X);
    P = E';
    
    % Check if the user specified the number of modes.
    if user_modes ~= 0
        num_modes = user_modes;
        total_var_z = 0;
    else 
        % If not specified, compute the principal components.
        z = (P * X)';
        
        % Compute the variance of each principal component.
        var_z = var(z, 1);
        
        % Calculate the total variance.
        sum_var_z = sum(var_z);
        
        % Compute the percentage of variance explained by each mode.
        total_var_z = (var_z ./ sum_var_z) * 100;
      
        % Retain the number of modes required to explain at least 95% of the total variance.
        for i = 1:10
            if sum(total_var_z(1:i)) >= 90
                break;
            end
        end
        num_modes = i;
    end

    % Retain only the first num_modes principal components.
    P = P(1:num_modes, :);
    z = (P * X)';
    
    % Normalize the retained principal components to be in the range [0, 1].
    x0 = 0; x1 = 1;
    for p = 1:num_modes
        z_ub = max(z(:, p));
        z_lb = min(z(:, p));
        z_norm = (x1 - x0) .* (z(:, p) - z_lb) ./ (z_ub - z_lb) + x0;
        temp(:, p) = z_norm';
    end
    
    z = temp;
end