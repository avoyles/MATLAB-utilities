function data_mat  = fitzpeaks_parser(filename_list, key_energies, gamma_lines, EoB_Time, attenuation_data, rhodrs, unc_rhodrs)
%Fitzpeaks_Parser   Reads in a list of Fitzpeaks reports, to parse and generate gamma-ray spectral analysis



% Interpolate attenuation coefficients
attenuation_model = pchip(attenuation_data(:,1).*1e3,attenuation_data(:,2));

% Set up cell array to hold results
data_mat = GammaCounts(length(key_energies),  length(filename_list));
% data_mat = cell([length(key_energies)  1+length(filename_list)]);
% data_mat(:,1) = num2cell(key_energies);

% Set up arrays to hold shelf and detector info
shelves = zeros(1,length(filename_list));
detectors = shelves;

% Update - now wrapping all counts into GammaCounts!  Construct as:
% obj = GammaCounts( E_gamma, t_half, number_of_counts, unc_number_of_counts, I_gamma, unc_I_gamma, live_time, mu, rhodr, file_name, mass, EoB_Time, unc_rhodr, covariance_data)


% Iterate over all Fitzpeaks files
for i=1:length(filename_list)
    
    % Open and read files
    fname = char(filename_list(i));
    fid = fopen(fname);
    
    % Get raw text for header regex extraction, find number of heder lines
    raw_str = fileread(fname);
    whichline = 1+find(contains(regexp(raw_str,'\n','split'),'Signif'));
    
    % Parse column data to cell structure
    %   'whichline' refers to the number of header lines in the report file
    parsed_fitzpeaks_data = textscan(fid,'%f %f %f %f %f %d %f %f %f %d','headerlines',whichline);
    fclose(fid);
    
    % Get shelf position, for efficiency correction
    shelf = regexp(raw_str, 'Shelf:\s+(\d+)[on]?', 'tokens');
    shelves(i) = cell2mat(cellfun(@(x) str2double(x{:}), shelf, 'UniformOutput', false));
    shelf = shelves(i);
    
    % Get detector ID, for efficiency correction
    detector = regexp(raw_str, 'Detector:\s+(\d+)[on]?', 'tokens');
    detector = cell2mat(cellfun(@(x) str2double(x{:}), detector, 'UniformOutput', false));
    if isempty(detector)
        detector = 1;
    end
    detectors(i) = detector;
    
    
    
    
    %     Select efficiency curve based on shelf
    if detector==1
        if shelf==5
            effcal = 'effcal_05.mat';
        elseif shelf==10
            effcal = 'effcal_10.mat';
        elseif shelf==15
            effcal = 'effcal_15.mat';
        elseif shelf==18
            effcal = 'effcal_18.mat';
        elseif shelf==22
            effcal = 'effcal_22.mat';
        elseif shelf==30
            effcal = 'effcal_30.mat';
        elseif shelf==40
            effcal = 'effcal_40.mat';
        elseif shelf==401
            effcal = 'effcal_40alt.mat';
        elseif shelf==50
            effcal = 'effcal_50.mat';
        elseif shelf==501
            effcal = 'effcal_50.mat';
        elseif shelf==70
            effcal = 'effcal_70.mat';
        elseif shelf==80
            effcal = 'effcal_80.mat';
        end
    elseif detector==2
        if shelf==1
            effcal = eff_1;
        end
    elseif detector==3
        if shelf==1
            effcal = eff_1;
        end
    elseif detector==4
        if shelf==1
            effcal = eff_1;
        end
    elseif detector==5
        if shelf==1
            effcal = eff_1;
        end
    elseif detector==6
        if shelf==1
            effcal = eff_1;
        end
    elseif detector==7
        if shelf==1
            effcal = eff_1;
        end
    end
    
    
    % Number of columns for temporary storage matrix
    num_cols = 9;
    
    % Get list of gamma energies from the current report
    extract_energies = parsed_fitzpeaks_data{1,1};
    dim = length(extract_energies);
    
    % Make matrix to hold data for this spectrum
    sto_mat = zeros(dim,num_cols);
    
    % Pull out gamma-ray energies from Fitzpeaks report
    sto_mat(:,1) = extract_energies;
    
    
    
    
    % Extract header
    mass = regexp(raw_str, 'Mass:\s+(\d+)', 'tokens');
    mass = cell2mat(cellfun(@(x) str2double(x{:}), mass, 'UniformOutput', false));
    mass_str = num2str(mass);
    foil_id = str2num(mass_str(1:end-2));
    
    ct = regexp(raw_str, 'Datetime:\s+([-\w: ]+)', 'tokens');
    ct = juliandate(datetime(ct{1,1}{1,1}));
    
    cl = regexp(raw_str, 'Live:\s+(\d+)', 'tokens');
    cl = cell2mat(cellfun(@(x) str2double(x{:}), cl, 'UniformOutput', false));
    
    % Append data from header to columns
    sto_mat(:,4) = mass;
    sto_mat(:,5) = ct;
    sto_mat(:,6) = cl;
    
    
    % Pull out gamma-ray energies and net counts
    %                     energy                     net counts                               %error
    sto_mat(:,1:3) = [extract_energies      double(parsed_fitzpeaks_data{1,6})      parsed_fitzpeaks_data{1,7}];
    
    
    % Replace energy values with closest from key energies, otherwise delete?
    % Indices of rows to keep
    ind = zeros(2,dim);
    for j=1:dim
        [Ydiff, idx] = min(abs(key_energies - sto_mat(j,1)));
        if Ydiff <= 1.0
            ind(:,j) = [j;idx];
        end
    end
    
    % Remove null values
    % ind = ind(ind ~= 0);
    ind = ind(:,ind(1,:) ~= 0);
    
    % Extract rows within 1 keV of key energies
    sto_mat = sto_mat(ind(1,:), :);
    
    % Replace with key energies
    sto_mat(:,1) = key_energies(ind(2,:));
    % append the t12 and Igamma data from a key
    sto_mat(:,7:9) = gamma_lines(ind(2,:),:);
    
    
    
    % Build up set of GammaCounts objects, using the parsed Fitzpeaks data.
    mu = ppval(attenuation_model,sto_mat(:,1));
    for gamma_index = 1:length(sto_mat(:,1))
        %  Construct as:
        % obj = GammaCounts( E_gamma, t_half, number_of_counts, unc_number_of_counts, I_gamma, unc_I_gamma, live_time, mu, rhodr, file_name, mass, EoB_Time, unc_rhodr, covariance_data)
        temp_object = GammaCounts( sto_mat(gamma_index,1), sto_mat(gamma_index,7), sto_mat(gamma_index,2), sto_mat(gamma_index,3), sto_mat(gamma_index,8), sto_mat(gamma_index,9), sto_mat(gamma_index,6), mu(gamma_index), rhodrs(foil_id), fname, sto_mat(gamma_index,4), sto_mat(gamma_index,5), EoB_Time, unc_rhodrs(foil_id), effcal);
        % Merge matrix into data_mat
        data_mat(ind(2,gamma_index),i) = temp_object;
%         data_mat(ind(2,gamma_index),i+1) = {temp_object};
    end
    
    
    
end

% Print unique shelf positions, for info
fprintf('Unique shelf positions:')
fprintf(' %2d',unique(shelves))
fprintf('\n')

fclose('all');

end

