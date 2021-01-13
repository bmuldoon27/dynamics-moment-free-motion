%%  Sandbox for Moment Free Motion of a Rigid Body
clc; clear all; close all;
%% Define Geometry Cases
% Symmetric Body
a1 = 2;
b1 = 2;
c1 = 2;
dims1 = [a1,b1,c1];

% Axisymmetric Body
a2 = 2;
b2 = 2;
c2 = 7;
dims2 = [a2,b2,c2];

% Asymmetric Body
a3 = 2;
b3 = 4;
c3 = 7;
dims3 = [a3,b3,c3];
%% Set Final Simulation Time
T = 1; % seconds
%% Define  Initial Conditions
% Initial mass center velocity:
xdot0 = [0, 0, 20];          % [x1dot(0), x2dot(0), x3dot(0)], m/s
% Initial mass center position:
x0 = [0, 0, 0];             % [x1(0), x2(0), x3(0)], m
% Initial Conditions for Quaternion Components (Orientation):
phi0=pi;
e00=cos(phi0/2); 
e10= 0.2;
e20= sqrt(1-e10^2);
e30= 0;
p0 = [e00, e10, e20, e30]; 
% Initial angular velocity      % [omega1(0), omega2(0), omega3(0)], m
omega0 = [10, 10, 10];
% Initial angular acceleration  % [omegadot1(0), omegadot2(0), omegadot3(0)], m
omegadot0 = [0, 0, 0];
% Compute initial rates of change of quaternions from prescribed omega0
% ---> Utilizing the quaternion representation of the angular velocity
L0 = [-e10, e00, e30, -e20;
      -e20, -e30, e00, e10;
      -e30, e20, -e10, e00];
pdot0 = 1/2*L0'*omega0';
% Construct the intial condition vector for numerical solver
IC = [xdot0, x0, omega0, p0, omegadot0, pdot0'];

%% Generate Solutions
[fig1_sym, fig2_sym, fig3_sym] = MMFMRB_quat(dims1,IC,T);
saveas(fig1_sym,'Sym_Sol_Params.png')
saveas(fig2_sym,'Sym_Sol_ICs.png')
saveas(fig3_sym,'Sym_Sol_Animation.png')
% 
% [fig1_axisym, fig2_axisym, fig3_axisym] = MMFMRB_quat(dims2,IC,T);
% saveas(fig1_axisym,'Axisym_Sol_Params.png')
% saveas(fig2_axisym,'Axisym_Sol_ICs.png')
% saveas(fig3_axisym,'Axisym_Sol_Animation.png')
% 
% [fig1_asym, fig2_asym, fig3_asym] = MMFMRB_quat(dims3,IC,T);
% saveas(fig1_asym,'Asym_Sol_Params.png')
% saveas(fig2_asym,'Asym_Sol_ICs.png')
% saveas(fig3_asym,'Asym_Sol_Animation.png')

