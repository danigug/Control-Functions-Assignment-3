clear all; close all; clc
%% Parameters
m = 1575; % [Kg] vehicle_mass
Jz = 2875; % [kg*m^2] yaw_mass_moment_inertia
a = 1.3; % [m] front_semi_wheelbase
b = 1.5; % [m] rear_semi_wheelbase
l = a+b;
CF = 2*60000; % [N/rad] front_axle_cornering_stiffness
CR = 2*60000; % [N/rad] rear_axle_cornering_stiffness

delta_max = 25; % [deg] maximum_steering_angle

Vv = 80/3.6; % [m/s]
%% State Space definition
% x = [beta, psi_dot]
% u = [delta_f, delta_r]
% y = [beta, psi_dot, rho_G, alfa_f, alfa_r, a_y]
A=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
    (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
B=[CF/(m*Vv) CR/(m*Vv);
    (CF*a/Jz) -(CR*b/Jz)];
C = [1,0
    0,1
    (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
    -1, -a/Vv
    -1, b/Vv
    (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
D = [0 0;
    0 0;
    CF/(m*Vv^2) CR/(m*Vv^2)
    1 0
    0 1
    CF/m CR/m];
% G = ss(A,B,C,D);
V=Vv;

A_c = [0 1 0 0
    0 -(CF+CR)/(m*V) (CF+CR)/m (CR*b-CF*a)/(m*V)
    0 0 0 1
    0 (CR*b-CF*a)/(Jz*V) (CF*a-CR*b)/Jz -(CR*b^2+CF*a^2)/(Jz*V)];

B1_c = [0
    CF/m
    0
    (CF*a)/Jz];

B2_c = [0
    ((CR*b-CF*a)/(m*V))-V
    0
    -(CR*b^2+CF*a^2)/(Jz*V)];

K=place(A_c,B1_c,[-1 -10 -50 -100])
Kff = ((m*V^2)/l)*(b/CF-a/CR+(a*K(3)/CR))+l-l*K(3)

sim("model.slx")

%% Results plots

% TO DO

