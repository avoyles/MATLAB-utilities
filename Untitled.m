clear all; clc;

t=linspace(0,30,1000);

y=-1.5+1.5.*cosh(1.11.*t)+7.38.*exp(-1.5.*t) .*sinh(1.11.*t);

plot(t,y)