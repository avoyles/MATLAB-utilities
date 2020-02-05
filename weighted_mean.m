function [mean, sigma_mean] = weighted_mean(x, sigma)
%weighted_mean  - Calculates weighted mean of a vector of values
%   Inputs:
%   x:     measurement values
%   sigma: uncertainty in each measurement value x_i
% 
%   Outputs:
%   mean:        Uncertainty-weighted mean
%   sigma_mean:  Propagated uncertainty in weighted mean
% 
% 
%   Assumes a vector of measurements:   x_i +/-  sigma_i
% 
%   Weights are calculated as:  w_i = 1 / (sigma_i)^2 
% 
% 
%
% See https://physics.stackexchange.com/a/329412
%     https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
% 
%
% A. Voyles
% 31 Oct 2017 - Version 0.1

% Catch N > lambdas
if length(sigma) ~= length(x)
    error('Vectors are of different lengths. \n %i measurement values provided. \n  %i uncertainties provided.', length(x), length(sigma))
end

w = 1./(sigma.^2);

mean = sum(x .* w) / sum(w);

sigma_mean = sqrt(1 / sum(w));


end

