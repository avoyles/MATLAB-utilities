function output = fraction_convert(fraction, masses, basis)

% Catch unequal # of species input
if length(fraction) ~= length(masses)
    error('Unequal number of species fractions and masses input. \n %i species fractions input. \n %i species molar masses input.', length(fraction), length(masses))
%     error(strcat(length(fraction),' number of fractions input.'))
%     error(strcat(length(masses),' number of molar masses input.'))
end

if strcmpi(basis, 'a')
    % Input is atomic percents
    output = (fraction .* masses) / sum(fraction .* masses);
    output = [output; output.*100];
elseif strcmpi(basis, 'w')
    % Input is weight percents
    output = (fraction ./ masses) / sum(fraction ./ masses);
    output = [output; output.*100];
else
    % Input is unknown - error out
    error('Unrecognized fraction basis flag')
end

end