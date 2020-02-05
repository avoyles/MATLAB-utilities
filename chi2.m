function output = chi2(y,y2,sigma)

% See http://maxwell.ucsc.edu/~drip/133/ch4.pdf

% % Catch unequal # of species input
% if length(fraction) ~= length(masses)
%     error('Unequal number of species fractions and masses input. \n %i species fractions input. \n %i species molar masses input.', length(fraction), length(masses))
% %     error(strcat(length(fraction),' number of fractions input.'))
% %     error(strcat(length(masses),' number of molar masses input.'))
% end

top = sum(y./sigma.^2);

bottom = sum(1./sigma.^2);

ybar = top./bottom;

output = sum( (y-y2).^2./sigma.^2   );

end