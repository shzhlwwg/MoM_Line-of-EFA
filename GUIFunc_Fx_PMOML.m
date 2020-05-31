% a = 100e-6;
% s = 100e-6;
% d = 100e-6;
% N = 100;
% M = 10;
% Ndots = 10;
% Coef_fric = 0.1;
function [MThetax,MPhi,MFx]=GUIFunc_Fx_PMOML(a,d,s,N,M,Ndots,Coef_fric)

p =a+s; % the pitch
epslon = 8.8514e-12;
k = 1/4/pi/epslon;


Thetax = linspace(-pi,pi,Ndots);
Phid = linspace(-pi,pi,Ndots);
Phi = -0.5*Thetax + Phid;
[MThetax,MPhi] = meshgrid(Thetax,Phi);

Fx = zeros(Ndots,Ndots);
Fy = zeros(Ndots,Ndots);
MFx = zeros(Ndots,Ndots);

for ind_thetax = 1:Ndots
    for ind_phi = 1:Ndots
        thetax = MThetax(ind_thetax,ind_phi);
        phi = MPhi(ind_thetax,ind_phi);
        
        dx = thetax/(2*pi)*3*p;
        % Input three-phase voltages
        % Vm = [0 1 0 0 1 0]';
        Vm = [sin(phi) sin(phi-2*pi/3) sin(phi+2*pi/3) sin(phi+2*pi/3) sin(phi-2*pi/3) sin(phi)]';
        
        %% Periodic MOM line Model
        %  the potential formula: u = 1/(2*pi*epslon)*log(r0/r)*q
        %                            y
        %                            |
        %            ---1---      ---2---      ---3---
        % ---------------------------|-----------------------> x
        %            ---4---      ---5---      ---6---
        %                            |
        R = a/N/2; % the radius of elements
        % r0(x0,y0), calculate the postion of line elements
        x0 = zeros(2*3*N,1);
        y0 = zeros(2*3*N,1);
        for i=1:2*3
            for j= 1:N
                y0((i-1)*N+j) = (3-i+0.5)/abs(3-i+0.5)*d/2; % y coordinates
                if i <= 3
                    x0((i-1)*N+j) = (j-(1+N)/2)*(2*R) + (i-2)*p + 1/2*dx;
                else
                    x0((i-1)*N+j) = (j-(1+N)/2)*(2*R) + (i-5)*p - 1/2*dx;
                end
            end
        end
        
        % r matrix: relative position matrix of every line to others.
        % rMatrix = zeros(2*3*N, 2*3*N);
        xMatrix = zeros(2*3*N, 2*3*N);
        yMatrix = zeros(2*3*N, 2*3*N);
        for i = 1:2*3*N
            for j = 1:2*3*N
                %         if i == j
                % %             rMatrix(i,j) = R;
                %         else
                xMatrix(i,j) = x0(i)-x0(j);
                yMatrix(i,j) = y0(i)-y0(j);
                %             rMatrix(i,j) = sqrt(xMatrix(i,j).^2 + yMatrix(i,j).^2);
                %         end
            end
        end
        % D matrix
        Dmatrix = zeros(2*3*N, 2*3*N);
        for i = 1:2*3*N
            for j = 1:2*3*N
                Dmatrix(i,j) =  Func_PMOML(xMatrix(i,j),yMatrix(i,j),R,M,p);
            end
        end
        
        %% A matrix
        A = zeros(2*3,2*3*N);
        for i = 1:2*3
            A(i,(i-1)*N + 1:i*N) = ones(1,N);
        end
        
        %% voltage
        Volt = A'*Vm;
        %% Q
        q = Dmatrix\Volt;
        q0 = q;
        
        % smooth
        for i = 1:2*3
            qa = q((i-1)*N+1:i*N);
            qq = q((i-1)*N+1:i*N);
            % qa(2);
            qq(2) = 1/2*qa(2)+1/4*(qa(1)+qa(3));
            %qa(4)
            qq(4) = 1/2*qa(4)+1/4*(qa(3)+qa(5));
            %qa(5)
            qq(5) = 1/2*qa(5)+1/4*(qa(4)+qa(6));
            % qa(n-1);
            qq(end-1) = 1/2*qa(end-1)+1/4*(qa(end)+qa(end-2));
            %qa(n-3)
            qq(end-3) = 1/2*qa(end-3)+1/4*(qa(end-2)+qa(end-4));
            %qa(n-4)
            qq(end-4) = 1/2*qa(end-4)+1/4*(qa(end-3)+qa(end-5));
            %qa(3)
            qq(3) = 1/16*(qa(1)+4*qa(2)+6*qa(3)+4*qa(4)+qa(5));
            %qa(n-2)
            qq(end-2) = 1/16*(qa(end)+4*qa(end-1)+6*qa(end-2)+4*qa(end-3)+qa(end-4));
            q((i-1)*N+1:i*N) = qq;
        end
        
        %% Coulombian force
        fx = zeros(3*N,1);
        fy = zeros(3*N,1);
        
        for i=1:3
            for j = 1:N
                % r1(x1,y1) relative postion between two charge
                x1 = zeros(3*N,1);
                y1 = zeros(3*N,1);
                fxx = zeros(3*N,1);
                fyy = zeros(3*N,1);
                for ii = 1:3
                    for jj = 1:N
                        
                        x1((ii-1)*N+jj) = x0((i-1)*N+j) - x0((ii-1+3)*N+jj);
                        y1((ii-1)*N+jj) = y0((i-1)*N+j) - y0((ii-1+3)*N+jj);
                        [fxx((ii-1)*N+jj),fyy((ii-1)*N+jj)] = Func_Exy_PMOML(x1((ii-1)*N+jj),y1((ii-1)*N+jj),p,M);
                        
                        fxx((ii-1)*N+jj) = fxx((ii-1)*N+jj)*q0((i-1)*N+j).*q0((ii-1+3)*N+jj);
                        fyy((ii-1)*N+jj) = fyy((ii-1)*N+jj)*q0((i-1)*N+j).*q0((ii-1+3)*N+jj);
                    end
                end
                
                fx((i-1)*N+j) = sum(fxx);
                fy((i-1)*N+j) = sum(fyy);
            end
        end
        
        Fx(ind_thetax,ind_phi) = sum(fx);
        Fy(ind_thetax,ind_phi) = sum(fy);
        
        if Fy(ind_thetax,ind_phi) >=0
            MFx(ind_thetax,ind_phi)= Fx(ind_thetax,ind_phi);
        else
            if Fx(ind_thetax,ind_phi)>0
                MFx(ind_thetax,ind_phi) = Fx(ind_thetax,ind_phi) + Coef_fric*Fy(ind_thetax,ind_phi);
            else
                MFx(ind_thetax,ind_phi) = Fx(ind_thetax,ind_phi) - Coef_fric*Fy(ind_thetax,ind_phi);
            end
        end
        
    end
end
MFx = MFx./3/p;


