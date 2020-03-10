function [tendl_E,tendl_xs] = tendl_plot(plotted_energies, tendl_path,abundances,varargin)
%tendl_plot Takes a list of isotopic abundances, to construct
%natural-abundnace target cross sections from monoisotopic TENDL output


% disp("Number of input arguments: " + nargin)
% disp("Number of optional input arguments: " + length(varargin))

% varargin{1}
% test_path

% Check to verify that number of isotopes = number of tendl files
if length(abundances) ~= length(varargin)
    fprintf('WARNING: Number of isotopic abundances (%d) does not match number of TENDL files (%d) supplied!\n',length(abundances), length(varargin))
    return
end

% Set energy grid for plotting
xp = plotted_energies;


% Get data structure from first file
tendl_data = dlmread(strcat(tendl_path,varargin{1}),' ',5,1);
tendl_xs = (abundances(1).*tendl_data(:,2));
tendl_E  = tendl_data(:,1);


if length(abundances) > 1
    % Loop over all remaining isotopes
    for i=2:length(abundances)
        tendl_data = dlmread(strcat(tendl_path,varargin{i}),' ',5,1);
        tendl_xs = tendl_xs + (abundances(i).*tendl_data(:,2));
    end
    
%     x = dlmread('./tendl/deuterons/235U/rp093236.L01.txt',' ',5,1);
%     y = dlmread('./tendl/deuterons/238U/rp093236.L01.txt',' ',5,1);
%     tendl_E  = x(:,1);
%     tendl_xs = (frac_natU235.*x(:,2)) + (frac_natU238.*y(:,2));
end


tendl_xs = spline(tendl_E, tendl_xs, xp);
tendl_E=xp;


end

