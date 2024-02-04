% The intgasp function models two nested functions, f and g.
function [model_f,model_g] = intgasp(design_f,response_f,design_g,response_g)

% design_f -- design and response for the "inside" function.
% design_g -- independent design and response for the "outside" function.
% response_f -- output of function f
% response_g -- output of funciton g
% The output of f (response_f) is used as an additional input to g.


% Initialize structures to hold the models for f and g.
model_f=struct();
model_g=struct();

% Set a default value for the hyperparameter alpha for both models. 
% It is used in the squared exponential kernel function.
model_g.alpha=2.0;
model_f.alpha=2.0;

% Assign input and output data for function f.
model_f.input=design_f;    % Independent inputs for f
model_f.output=response_f; % Outputs of f, which are also used as dependent inputs to g.

% Assign input and output data for function g.
model_g.input=design_g;    % Independent inputs for g.
model_g.output=response_g; % Outputs of g.

% Dimensionality of the input and output spaces for f and g.
dimu=size(model_f.input,2); % Dimension of the inputs to f.
dimy=size(model_g.input,2); % Dimension of the independent inputs to g.
dimz=size(model_f.output,2);% Dimension of the outputs of f (dependent inputs to g).

% Check for a size constraint: the number of output points from f should match 
% the number of input points to g.
if(size(model_g.output,1) ~= size(model_f.input,1))
    error('the number of rows in the input should be the same as the number of the rows in the output');
end

% Configuration for modeling function f with ppgasp().
% Set the kernel type, and determine the dimensions for kernel and hyperparameter configurations.
k_type_f=["pow_exp"];
dim_io_f = size([model_f.input],2);
optionsf.kernel_type=repmat(k_type_f,dim_io_f,1);
optionsf.alpha=repmat(model_f.alpha,dim_io_f,1);
optionsf.trend=[ones(size(model_f.input,1),1) model_f.input];
optionsf.zero_mean = false;
optionsf.nugget_est = true;

% GP regression of  f using ppgasp().
model_f = ppgasp(model_f.input, model_f.output, optionsf);
model_f.tau = model_f.nugget; % Store the estimated nugget as tau.

% Configuration for modeling function g with ppgasp().
% Since g takes inputs from both its own independent input set and the outputs of f, 
% it has a (potentially) larger input space.
k_type=["pow_exp"];
dim_io=size([model_g.input,model_f.output],2);
options.kernel_type=repmat(k_type,dim_io,1);
options.alpha=repmat(model_g.alpha,dim_io,1);
options.trend=[ones(size(model_g.input,1),1)  [model_g.input, model_f.output]];
options.zero_mean = false;
options.nugget_est = true;

% Use Gaussian Process Regression to model function g.
model_g=ppgasp([model_g.input, model_f.output], model_g.output, options);

end