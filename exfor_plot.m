function h = exfor_plot(matrix,varargin)
%weighted_mean  - Calculates weighted mean of a vector of values
%   Inputs:
%   matrix:     4-column matrix of EXFOR data
%               X         +-dX         Y	        +-dY
%               MeV       MeV       barns	       barns
% 
%   varagin:  optional arguments for plot options
%
%
% 
%   Outputs:
%   h:        handle for plot
% 
% 
%
% A. Voyles
% 07 Nov 2017 - Version 0.1

% Catch wrong matrix dimensions
[~,n] = size(matrix);

if n ~= 4
    error('Matrix is of wrong dimension. \n Matrix has %i columns of data, and 4 are expected.', n)
end

% Check for zero uncertainties
if all(matrix(:,4) ==0)
    % All uncertainties are equal to zero
    h=plot(matrix(:,1),1e3.*matrix(:,3),varargin{:});
else
    % Actual uncertainties reported
    h=errorbar(matrix(:,1),1e3.*matrix(:,3),1e3.*matrix(:,4),1e3.*matrix(:,4),matrix(:,2),matrix(:,2),varargin{:});
end

end

