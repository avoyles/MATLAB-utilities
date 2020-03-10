function [exfor_data,exfor_dataset_indices, exfor_names] = exfor_read(filename)
%exfor_read Parser to read EXFOR files and extract individual data sets for
%plotting

fid = fopen(filename);

raw_str = fileread(filename);
whichline = find(contains(regexp(raw_str,'\n','split'),'EXFOR-ID'));

% Parse column data to cell structure
%   'whichline' refers to the number of header lines in the report file
parsed_exfor = textscan(fid,'%f %f %f %f %s %s %s %s %s','headerlines',whichline);
[number_of_parsed_lines , ~] = size(parsed_exfor{1,1});
if number_of_parsed_lines==0
    fprintf('WARNING: Unable to read in file %s\n',fname)
end


unique_datasets = unique(parsed_exfor{1, 8}  );
number_of_unique_datasets = length(unique_datasets);
exfor_dataset_indices = zeros(2,number_of_unique_datasets);



for i=1:number_of_unique_datasets
    idx = find(strcmp([parsed_exfor{1, 8}], unique_datasets{i})); % single line engine
    exfor_dataset_indices(1,i) = min(idx);
    exfor_dataset_indices(2,i) = max(idx);
    exfor_dataset_indices;
end

exfor_names = parsed_exfor{1, 6}(exfor_dataset_indices(1,:));

for j=1:length(exfor_names)
    temp_string = exfor_names{j};
    year_string = temp_string(1:4);
    temp_string = erase(temp_string,',');
    temp_string = erase(temp_string,'+');
    temp_string = erase(temp_string,year_string);
    temp_string = strcat(temp_string,' (',year_string,')');
    
    exfor_names{j} = temp_string;
end

exfor_data = [parsed_exfor{1, 1:4}];

fclose(fid);



end

