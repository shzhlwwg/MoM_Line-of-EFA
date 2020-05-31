clear all;
close all;
% 最终程序, 调整 an,sn,dn 范围可以生成不同参数范围的电机输出力的数据
% Calculate the max/min Fx, Fy of EFAs at different a, s, d parameters,
% without consideration of breakdown
% Example parameter
a = 100e-6;
d = 100e-6;
s = 100e-6; %space

N = 100; % even, element number
M = 100;% array number
%
Na = 31;
Ns = 31;
Nd = 7;
an = [1e-6,linspace(10e-6,300e-6,Na-1)];
sn = [1e-6,linspace(10e-6,300e-6,Ns-1)];
dn = [10e-6,linspace(50e-6,300e-6,Nd-1)];

[An, Sn, Dn] = meshgrid(an,sn,dn);
Xmax = zeros(Ns,Na,Nd);
Ymax = zeros(Ns,Na,Nd);
Ymin = zeros(Ns,Na,Nd);

for k = 1:Nd
    for m = 1:Na
        for n = 1:Ns
            tic % tic 
            disp([n,m,k]) % disp process
            
            a = An(n,m,k);
            s = Sn(n,m,k);
            d = Dn(n,m,k);

            p=a+s;
            
            Num = 11;
            % thetax:0, pi/3:pi/2
            % phi : pi/4, -pi/3 pi/6
            % the most possible points where max or min forces generate
            thetax = [0 0 pi/3:pi/48:pi/2]'; 
            phi = [-pi/3 pi/6, pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 pi/4 ]';
            
            % Thetax = 2*pi/3/p*Dx;
            dx = 3*p/2/pi*thetax;
            
            Fcoulombian = zeros(Num,2);
            Q = zeros(Num,6);
            parfor i = 1:Num % parallel calculate to accelerate, needless if only one point.
                [Fcoulombian(i,:),Q(i,:)] = Func_MatrixChargeCoulF_PMOML(a,d,s,dx(i),phi(i),N,M);
            end
            % force per area per voltage
            Fcx = Fcoulombian(:,1)/3/p;
            Fcy = Fcoulombian(:,2)/3/p;
            Xmax(n,m,k) = max(Fcx);
            Ymax(n,m,k) = max(Fcy);
            Ymin(n,m,k) = min(Fcy);
            
            disp([toc])    % disp time costed every point
        end
    end
end
%% save mat data
% save('Coul_PMOML_A_S_D50-100.mat','An','Sn','Dn','Xmax','Ymax','Ymin')
% surfc(An(:,:,2),Sn(:,:,2),Xmax(:,:,2))
% xlim([0,3]*1e-4)
% ylim([0,3]*1e-4)
%% write xls data
% for i = 1:7
%     A = An(:,:,i);
%     S = Sn(:,:,i);
%     Xa = Xmax(:,:,i);
%     Ya = Ymax(:,:,i);
%     Yi = Ymin(:,:,i);
%     
%     xlswrite('Coul_PMOML_A1-300_S1-300_D10-300.xlsx',...
%         [A(:),S(:),Xa(:),Ya(:),Yi(:)],num2str(dn(i)*1e6));
% end

