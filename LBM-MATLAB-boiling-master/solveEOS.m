%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solveEOS.m: to solve the EOS for the coexisting densities 
%                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Shan Chen Lattice Boltzmann sample in Matlab
% Reference: Li, Qing, et al. "Lattice Boltzmann modeling of boiling heat 
%            transfer: The boiling curve and the effects of wettability." 
%            International Journal of Heat and Mass Transfer 85 (2015): 
%            787-796.

clear, clc;

a   = 3/49;
b   = 2/21;
R   = 1;
ome = 0.344;
Tc  = 0.0778/0.45724*a/(b*R);
T   = 0.86*Tc;
V   = [0.15:0.0001:3.5];
phi = (1 + (0.37464 + 1.54226*ome - 0.26992*ome^2)*(1 - sqrt(T/Tc)))^2;
p   = R*T./(V - b) - a*phi./(V.^2 + 2*b.*V - b^2);

plot(V, p, 'b');
hold on;

title('Peng-Robinson EOS');
xlabel('Vm');
ylabel('P');

P = 0.025;
deltaP = 0.025/2;
S1 = 1;
S2 = 0;
S = (S1 - S2)/S1;
%while(abs(S) > 0.00000001)
while(abs(S)>0.0001)
    Vm = solve('R*T/(Vm - b) - a*phi/(Vm^2 + 2*b*Vm - b^2) - P', 'Vm');
    % to get solutions
    Vm1 = double(vpa(real(subs(Vm(1)))));
    Vm2 = double(vpa(real(subs(Vm(2)))));
    Vm3 = double(vpa(real(subs(Vm(3)))));
    V1  = min([Vm1 Vm2 Vm3]);
    V2  = max([Vm1 Vm2 Vm3]);
    
    fcs = @(Vm) R*T./(Vm - b) - a*phi./(Vm.^2 + 2*b.*Vm - b^2);
    S1  = integral(fcs, V1, V2);
    % S1 is the integration area
    S2 = P*(V2 - V1);
    % S2 is a rectangular area
    if(S1 > S2)
        P = P + deltaP;
    else
        P = P - deltaP;
    end
    deltaP = deltaP/2;
    P
    S = (S1 - S2)/S1;
    S
end
P
V1
V2
rho1 = 1/V2
rho2 = 1/V1
plot(V, P*ones(1, length(V)), 'r');
plot([V1 V2], [P P], 'ko');
hold off;
    