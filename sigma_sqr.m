function [sigma_sqr] = sigma_sqr(H,invR,y,theta, m,q) 
  % Computes the squared difference between observed output and predicted output.
  % This is the residual sum of squares (RSS).
  squared_diff = transpose(y - H*theta) * invR * (y - H*theta);
  
  % Normalizes the RSS by the difference between the number of observations
  % and the number of parameters (degrees of freedom).
  sigma_sqr = squared_diff ./ abs(m - q);
end