% kernel_PE function computes the kernel evaluations between two sets of input samples 
% using a power exponential kernel with specified alpha and parameters.
function R = kernel_PE(a,b,alpha,param)

    % Determine the number of samples in input sets a and b.
    n = size(a,1); % Number of samples in input a.
    m = size(b,1); % Number of samples in input b.
    
    % Determine the dimensionality of the inputs.
    A = size(a);
    d = A(2); % Dimension of the inputs.
    
    % Initialize the exponent matrix with zeros.
    exponent = zeros(n,m);
    
    % Loop through all combinations of samples from a and b.
    for k = 1:n
        for l = 1:m
            % Compute the distance between the k-th sample from a and the l-th sample from b 
            % in each dimension using the power-exponential kernel formula.
            for j = 1:d
                r(j) = (abs(a(k,j) - b(l,j)).^alpha)./(param(j).^2.0);
            end
            % Sum the computed distances for all dimensions.
            exponent(k,l) = -sum(r);
        end
    end
    
    % Apply the exponential function to the computed exponents to get the kernel evaluations.
    R = exp(exponent);

end