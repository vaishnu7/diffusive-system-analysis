clc;
clear;
close all;
format long
syms u v
a = 4;
b = 2;
d = 1.85;
h = 0.118;

jac = jacobian([u*(1-u) - a*u*v / (u + v) - h,v*(-d + b*u/u+v)], [u,v]); %we calculate Jacobian for non-diffusive system only

[uu, vv] = vpasolve( [u*(1-u) - a*u*v / (u + v) - h == 0, v*(-d + b*u/u+v) == 0], [u,v] );
pp = [uu, vv]
