function phi = rk4(fname,phi0,t)
% RUNGE-KUTTA 4TH ORDER INITAL VALUE ODE SOLVER
% Solves via the 4th-order Runge Kutta method the first-order ODE given by:
%
%   dphi   
%   ---- = f(phi,t)  with Initial Condition phi(t0) = phi0
%    dt
%
% INPUTS:
%       Right-Hand-Side Function: fname
%       Initial Point: phi0
%       Time Vector: t
% OUTPUTS:
%       Phi Column Vector: phi
% FORM:
%       phi = rk4(fname,phi0,t)
% Mike Hansen
% 5-2011

N = length(t);
Neq = length(phi0);
dt = t(2) - t(1);
phi = zeros(Neq,N);
phi(:,1) = phi0;
for i = 2:N
    dy1 = dt*feval(fname,phi(:,i-1),t(i-1));
    dy2 = dt*feval(fname,phi(:,i-1)+dy1,t(i-1)+dt/2);
    dy3 = dt*feval(fname,phi(:,i-1)+dy2,t(i-1)+dt/2);
    dy4 = dt*feval(fname,phi(:,i-1)+dy3,t(i-1)+dt);
    phi(:,i) = phi(:,i-1) + 1/6*(dy1 + 2*dy2 + 2*dy3 + dy4);
end