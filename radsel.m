function radius = radsel(X, P, D, N, STD, IQR)
%RADSEL Selects an optimal radius parameter for nonlinear measure estimation.
%   This function systematically selects an optimal radius parameter for 
%   estimating nonlinear measures derived from the correlation integral, 
%   utilizing the Kernel Density Estimator (KDE) framework. It is designed 
%   to enhance the accuracy of measures such as correlation dimension and 
%   Kolmogorov-Sinai entropy. The methodology is validated through numerical 
%   experiments on signals from nonlinear systems and experimental 
%   electroencephalographic time series.
%
%   RADIUS = RADSEL(X) returns the optimal radius for the input data matrix X, 
%   where rows represent samples and columns represent dimensions. The function 
%   uses default values for optional parameters.
%
%   RADIUS = RADSEL(X, P) specifies the order of the norm used in the calculation. 
%   Default is 2 (Euclidean distance).
%
%   RADIUS = RADSEL(X, P, D) considers D dimensions from the input data. The 
%   default is the number of columns in X.
%
%   RADIUS = RADSEL(X, P, D, N) considers N samples from the input data. The 
%   default is the number of rows in X.
%
%   RADIUS = RADSEL(X, P, D, N, STD) uses the standard deviation STD of the input 
%   data in optimizing the radius parameter. The default is the standard deviation 
%   computed across all dimensions of X.
%
%   RADIUS = RADSEL(X, P, D, N, STD, IQR) uses the interquartile range IQR of the 
%   input data. The default is the interquartile range computed across all dimensions 
%   of X.
%
%   Example:
%       % Load example data
%       load('sampleData.mat');
%       % Estimate optimal radius
%       optimalRadius = RADSEL(data);
%       disp(['Optimal Radius: ', num2str(optimalRadius)]);
%
%
%   References:
%   Johan Medrano, Abderrahmane Kheddar, Annick Lesne, Sofiane Ramdani; 
%   "Radius selection using kernel density estimation for the computation of 
%   nonlinear measures." Chaos, 1 August 2021; 31 (8): 083131. 
%   https://doi.org/10.1063/5.0055797.
%
%   Please cite the above work when using the RADSEL function for academic or 
%   research purposes.



% Setup arguments 
% =========================================================================

% Make X a column vector
% -------------------------------------------------------------------------
if size(X,1) == 1
    x = X'; 
else 
    x = X; 
end

% Order of the norm  
% -------------------------------------------------------------------------
try
    p = P; 
catch
    p = 2; 
end

% Number of dimensions 
% -------------------------------------------------------------------------
try
    d = D;
catch
    d = size(x, 2);
end

% Number of samples 
% -------------------------------------------------------------------------
try 
    n = N;
catch 
    n = size(x, 1);
end

% Standard deviation  
% -------------------------------------------------------------------------
try 
    s = STD; 
catch 
    s = std(x); 
end

% Interquartile range
% -------------------------------------------------------------------------
try 
    q = IQR; 
catch 
    q = iqr(x);
end 

% Calculate the different terms that compose the optimal radius
% =========================================================================

% Compute the scale
% -------------------------------------------------------------------------
scale = min(q / 1.34, s);

% Compute the density factor g
% -------------------------------------------------------------------------
g = n ^ (-1./(d + 4.)); 

% Volume of a ball in dimension d under norm p
% -------------------------------------------------------------------------
if p == 1 
    nballvolume = 2.^d / factorial(d); 
elseif p == 2
    nballvolume = sqrt(pi)^d / gamma(d / 2 + 1.); 
elseif p == inf
    nballvolume = 2.^d; 
elseif p > 0
    nballvolume = (2. * gamma(1 / p + 1.))^d / gamma(d / p + 1.); 
end

% Compute the ball correction factor a
% -------------------------------------------------------------------------
if d == 1
    a = (12 * sqrt(pi)) ^ (1./5.);
else 
    a = (( ... 
         4. * (2. * sqrt(pi))^d  * ...  
        (3. * gamma(1. + (d+2)/p) * gamma(1. + 1 / p))^2 ...
        ) / ( ... 
        (d + 2.) * nballvolume * ...
        (gamma(1 + 3./p) * gamma(1. + d/p))^2 ...
    )) ^ (1. / (d + 4.));
end

% Compute the optimal radius
% =========================================================================
radius = a * g * scale; 

end

