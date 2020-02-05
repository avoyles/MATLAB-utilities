function final_activity = bateman(N_0, N, lambdas, t_final)
%bateman  - Numerical solution to the Bateman equation
%   Returns final activity A_N(t) of Bateman Equation, for a simplified decay chain, at some time t=t_final:
%   1 -->  2 --> 3 --> 4 --> i --> N
%   N_1(t=0) = N_0
%   N_2(t=0) = N_3(t=0) = ... = N_N(t=0)  = 0
%
%   lambdas = [lambda_1  lambda_2  lambda_3  ...  lambda_N]
%
% See https://en.wikipedia.org/wiki/Radioactive_decay#Chain-decay_processes
% 
%
% A. Voyles
% 27 Feb 2016 - Version 0.1

% Catch N > lambdas
if length(lambdas) < N
    error('Insufficient decay constants provided. \n %i decay constants provided. \n Activity of chain index %i requested.', length(lambdas), N)
end

% Initialize activity
final_activity = 0;

%  Bateman solution:
%  A_N(t) = N_0 *  sum( c_i *exp(-lambda_i *t)   )

% Calculate c_i

% Initialize c_i,A_i
c_i = zeros(1,N);
A_i = c_i;

for j=1:N
    top = prod(lambdas(1:N));
    bottom = prod(nonzeros(lambdas(1:N)-lambdas(j)));
    % Catch bottom=0
    if bottom==0
        bottom=1;
    end
    c_i(j) = top./ bottom;
end


% Return final activity
for j=1:N
    A_i(j) = c_i(j) .*exp(-lambdas(j) .*t_final);
end
final_activity  = N_0 .*  sum( A_i   );

end

