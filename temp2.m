clear all; clc;

t=0:0.01:20;

v1i = 25;
v2i = 35;

x1i=45;
x2i=0;

a1=-2;
a2=-4;

disc = (v1i-v2i)^2 - .5*4*(x1i-x2i)*(a1-a2)

x1 = x1i + v1i.*t + .5.*a1.*t.^2;
x2 = x2i + v2i.*t + .5.*a2.*t.^2;

plot(t,x1,t,x2)
legend('Car1','Car2','location','best')