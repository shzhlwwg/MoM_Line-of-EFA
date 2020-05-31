function D = Func_PMOML(x,y,R,M,p)
% V = 2*k*SUM (log(r10/r20))*q
% periodic method of moments of line
% M must be even

epslon = 8.8514e-12;
k = 1/4/pi/epslon;

Ind = -M:M;
Dx = Ind*3*p;

r10 = sqrt(0^2+(0-Dx).^2);
r10(M+1) = R;
r20 = sqrt(y.^2+(x-Dx).^2);
r20(find(r20==0)) = R;
v = 2*k*log(r10./r20);

D = sum(v);


%%
%               |y
%               |
%               |
% -----o-----o--x--o-----o-----o----->x
%               |
%               |

% M = 100;
% epslon = 8.8514e-12;
% k = 1/4/pi/epslon;
% p = 200e-6;
% R = 10e-6;% the radius of lines
% % (xo,yo) the position of lines in global coordinate system
% xo = 0;
% yo = 0;
% % (x01,y01) zeros position
% x01 = 0;
% y01 = 0;
% % aim position (x,y)
% x = p;
% y = 0;
% 
% Ind = -(M-1):2:(M-1); % number :2*M
% Dx = Ind/2*3*p;
% % (xl,yl) positions of lines
% xl = Dx+xo;
% yl = yo;
% % r01 dist between zero point and lines
% r01 = sqrt((x01 - xl).^2+(y01-yl).^2);
% %r02 dist between point(X,Y) and the lines
% r02 = sqrt((y-yl).^2+(x-xl).^2);
% %(x,y) is too close to lines
% for m = 1:M
%     if r02(m) < R
%         r02(m) = R;
%     end
% end
% v = 2*k*log(r01./r02);
% D = sum(v)
