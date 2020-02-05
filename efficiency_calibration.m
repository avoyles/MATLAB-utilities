function [efficiency,unc_efficiency,pcov,popt] = efficiency_calibration(E_gamma,filename)
%effcal Generates an HPGE efficiency curve and 1-sigma uncertainty band
%   Detailed explanation goes here
%   See
%   https://github.com/amlewis1578/BernsteinGroupMeetings/blob/master/Slides/ALewis_Dec102018.pptx
%   for model formalism
% 
%   E_gamma:  gamma-ray energy
%   
%   popt:  5-parameters for efficiency curve
% 
%   pcov:  covariance matrix for the best-fit parameters
% 
% 

% Load popt and pcov from NPAT
fit_parameters=load(filename);
StrName=fieldnames(fit_parameters);
% Figure out which is which - pcov is n x n, popt is 1 x 5
[m, n] = size(fit_parameters.(StrName{1}));
if m==n
%     Then StrName{1} is pcov
    pcov = fit_parameters.(StrName{1});
    popt = fit_parameters.(StrName{2});
else
%     Then StrName{2} is pcov
    pcov = fit_parameters.(StrName{2});
    popt = fit_parameters.(StrName{1});
end



% beta = popt;
dim = length(popt);
num = length(E_gamma);

unc_efficiency = zeros(size(E_gamma));

% Define model formalism
eff_5_parameter = @(x,beta)(beta(1).*exp(-beta(2) .* (x.^beta(3)  )) .* (1 - exp(-beta(4) .* (x.^beta(5)  )))  );



efficiency = eff_5_parameter(E_gamma,popt);

% Numerically calculate the Jacobian matrix
% my_jacobian = zeros(length(E_gamma),dim);
my_jacobian = zeros(1,dim);
for j=1:num
    for i=1:dim
        % Define dBeta
        dBeta_i = popt(i) * 1E-8;
        
        % Create altered beta arrays
        tmp_popt_plus = popt; tmp_popt_plus(i) = popt(i) + (dBeta_i./2);
        tmp_popt_minus = popt; tmp_popt_minus(i) = popt(i) - (dBeta_i./2);
        
        % Use central difference approximation to df/dB_i
        my_jacobian(i) = ((eff_5_parameter(E_gamma(j),tmp_popt_plus) - eff_5_parameter(E_gamma(j),tmp_popt_minus)) ./ dBeta_i)';
    end
    
    % Uncertainty propagation
    unc_efficiency(j) = sqrt(my_jacobian * pcov * my_jacobian');
end



end

