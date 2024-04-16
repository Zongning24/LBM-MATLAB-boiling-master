%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization.m: initialize wall nodes and densities
%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shan Chen Lattice Boltzmann sample in Matlab
% Copyright Wei Gong
% Address: Nottingham NG7 2RD, UK
% E-mail: wei.gong@nottingham.ac.uk
% Reference: Li, Qing, et al. "Lattice Boltzmann modeling of boiling heat 
%            transfer: The boiling curve and the effects of wettability." 
%            International Journal of Heat and Mass Transfer 85 (2015): 
%            787-796.

global obst obst_b obst_u;
 
[y, x]    = meshgrid(1:ly, 1:lx); % lx=200; ly=150

% wall nodes and fluid nodes 壁结点和液结点
obst_b = (y == 1);                   %  obst=200*150
obst_u = (y == ly);                  %  obst=200*150
obst   = ((y == 1)|(y == ly));       %  obst=200*150   "|" 或者 
flui   = (ismember(y, 2:(ly - 1)));  %  flui=1*200     ismember 判断数组元素是否为集数组成员

% liquid-vopor system   液汽系统
rho = reshape(rho, lx, ly); % reshape 重构数组
rho(:, 1:88)   = rho_l; % 液体密度
rho(:, 92:ly)  = rho_v; % 气体密度
rho(:, 89:91)  = ones(lx, 1)*linspace(rho_l, rho_v, 3);

% % initial bubble 初始气泡
% rho(((x - lx/2).^2 + (y - ly/3).^2) < 20^2) = rho_v;
% rho(((x - lx/2).^2 + (y - 2*ly/3).^2) < 20^2) = rho_v;
% rho(((x - lx/4).^2 + (y - ly/3).^2) < 20^2) = rho_v;
% rho(((x - 3*lx/4).^2 + (y - ly/3).^2) < 20^2) = rho_v;
% rho(((x - lx/2).^2 + (y - 5*ly/100).^2) < 5^2) = rho_v;

rho = reshape(rho, 1, lxy);      % lxy=200*150
T   = reshape(T, lx, ly);
T(obst_u)      = Ts;             % 最顶层温度  Ts=0.86*Tc
T(flui)        = Ts;    
% T(obst_b)      = Ts;           % for debug use without temperature field
T(obst_b)      = Tb;             % Tb  = Ts + dT
% % small temperature fluctuations  小范围温度波动
T(1:4:lx, 2)   = T(1:4:lx, 2) + (linspace(0.3, 0, length(1:4:lx)))'*Ts;  
T(2:4:lx, 2)   = T(2:4:lx, 2) + (linspace(0.3, 0, length(2:4:lx)))'*Ts;
T(3:4:lx, 2)   = T(3:4:lx, 2) - (linspace(0.3, 0, length(3:4:lx)))'*Ts; 
T(4:4:lx, 2)   = T(4:4:lx, 2) - (linspace(0.3, 0, length(4:4:lx)))'*Ts;
% T((round(lx/2) - 2):(round(lx/2) + 2), 1) = Tb;
T = reshape(T, 1, lxy);
A([8 9], (rho > rho_a))  = 1/tau_l;
A([8 9], (rho <= rho_a)) = 1/tau_v;
lambda = 0.028*reshape(rho, lx, ly);
phi    = (1 + (0.37464 + 1.54226*ome - 0.26992*ome^2)*(1 - sqrt(Ts/Tc))).^2;
p      = rho.*R.*Ts./(1 - b.*rho) - ...
       a.*phi.*rho.^2./(1 + 2.*b.*rho - b^2.*rho.^2);

% distribution function 公式（5）
u_squ = ux.^2 + uy.^2;
me(1, :) = ones(1, lxy);
me(2, :) = -2 + 3*u_squ;
me(3, :) = 1 - 3*u_squ;
me(4, :) = ux; 
me(5, :) = -ux;
me(6, :) = uy;
me(7, :) = -uy;
me(8, :) = ux.^2 - uy.^2;
me(9, :) = ux.*uy;
me       = ones(9, 1)*rho.*me;

fe = reshape(Minv*me, 9, lx, ly);
ff = fe;

clear flui dT