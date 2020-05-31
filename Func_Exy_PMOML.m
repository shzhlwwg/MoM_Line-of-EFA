function [Ex,Ey] = Func_Exy_PMOML(x,y,p,NL)

% º∆À„÷‹∆⁄MoM
% a infinite line of charge: E = 2*k*q/r;
% a = 1;
% s = 1;
% NL = 100;
% 
% p = a+s;

epslon = 8.8514e-12;
k = 1/4/pi/epslon;
Ind = -NL:NL;
Dx = Ind*3*p;


Ex = sum(2.*k.*((x-Dx)./(y.^2+(x-Dx).^2)));
Ey = sum(2.*k.*((y)./(y.^2+(x-Dx).^2)));




