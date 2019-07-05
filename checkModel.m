close all
clear all

%%%%%%%%%%%%%% Check  Noplasma model CREATE for ISTTOK %%%%%%%%%%%%%%%%5


%%% Plasma-Less shot with heavisides in primary,vertical and horizontal
load('shot_45411.mat');
load('ISSTOK_nopl_newmesh_CL.mat');
%load('ISSTOK_nopl_CL.mat');

% Experimental
time = data.time*1e-6;
inputs = [data.prim, data.hor, data.vert];
Fluxes_real_1 = data.mirnv_corr/(50*49e-6)*1;
Fluxes_real = flip(Fluxes_real_1);

% % Equilibrium IOs
% yeq = LinearModel.YEquil(55:66)';
% ueq = LinearModel.XEquil(1:3)';

C = LinearModel.C(55:66,1:3);
% C = C([7:12 1:6],:);
% Fluxes = C*(inputs-ueq)' + yeq';
Fluxes = C*(inputs)';


% Eddy currents dynamics 
L11 = LinearModel.L(4:51,4:51);
L12 = LinearModel.L(4:51,1:3);
C1  = LinearModel.C([55:66],4:51);
C2  = LinearModel.C([55:66], 1:3); 
R   = LinearModel.R(4:51,4:51);
A   = -R*inv(L11);
B   = R*inv(L11)*L12;
C   = C1*inv(L11);
D   = C2-C1*inv(L11)*L12;
% C = C([7:12 1:6],:);
% D = D([7:12 1:6],:);
sys = ss(A,B,C,D);

y=lsim(sys,inputs,time);


%% Plot
figure(1)
 
subplot(1,2,2)
plot(inputs)
grid on
title('PFC currents')
legend('primary','Horizontal','Vertical')
 
 for i = 1 : 12
    subplot(1,2,1)
    cla
    plot(Fluxes(i,:))
    grid on
    hold on
    plot(Fluxes_real(i,:))
    hold on
    plot(y(:,i))
    legend('C*inputs','experimental','CREATE')
    title(['Probe #',num2str(i)])
    pause
 end
 
 
 
 
 
 
 
 %%%%
% A=LinearModel.A(1:51,1:51);
% A = -inv(LinearModel.L(1:51,1:51))*LinearModel.R(1:51,1:51);
% A(1:3,1:3)=-eye(3)*1e5;
% A(1:3,4:end)=0;
% B = inv(LinearModel.L(1:51,1:51))*LinearModel.S(1:51,1:end);
% B = B(:,1:3);
% B(1:3,1:3)=eye(3)*1e5;
% % B=LinearModel.B(1:51,1:3);
% C=LinearModel.C([66 78],1:51);
% theta = -15*pi/180;
% C=C(1,:)*cos(theta)+C(2,:)*sin(theta);
% D=LinearModel.D([66 78],1:3);
% D=D(1,:)*cos(theta)+D(2,:)*sin(theta);
% 
% C = [LinearModel.C(1:3,1:51); C];
% D = [LinearModel.D(1:3,1:3); D];

% %%%%%%
% pippo = ss(A,B,C,D);
% time = data.time*1e-6;
% y=lsim(pippo,inputs,time);

%%%%%%%
%%%%%%%