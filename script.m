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

Tsim = 1500;
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
Kff = ((m*V^2)/l)*(b/CF-a/CR+(a*K(3)/CR))+l-b*K(3)



%sim("model.slx")

%% Look-up table of control gains

B1_c = [0
    CF/m
    0
    (CF*a)/Jz];

AM = cell(1,13);
K_lookup = zeros(13,4);
Kff_lookup = zeros(13,1);

for i = 1:13 
    V(i) = (10*i)/3.6; % speeds from 10 to 130 Km/h expressed in m/s

    AM{i} = [0 1 0 0
    0 -(CF+CR)/(m*V(i)) (CF+CR)/m (CR*b-CF*a)/(m*V(i))
    0 0 0 1
    0 (CR*b-CF*a)/(Jz*V(i)) (CF*a-CR*b)/Jz -(CR*b^2+CF*a^2)/(Jz*V(i))];

    K_lookup(i,:) = place(AM{i},B1_c,[-1 -10 -50 -100]);
    Kff_lookup(i) = ((m*V(i)^2)/l)*(b/CF-a/CR+(a*K_lookup(i,3)/CR))+l-b*K_lookup(i,3);

end

%% Test

V = 55/3.6;

% sim("model.slx");

%% Evaluation along the relevant curvature profile
selector = 2;
close all

% TO DO: evaluate for different tunings

save_files = false;

V = 60/3.6;
Tsim = 30;

sim("model.slx");

name_fig = sprintf('Relevant curvature profile');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[s]'); ylabel('[m^-1]');
plot(tout,Kl_profile, 'LineWidth', 1)
axis tight;  
currentYlim = ylim;  
padding = 0.1 * (currentYlim(2) - currentYlim(1));  
ylim([currentYlim(1) - padding, currentYlim(2) + padding]);  % Set new y-axis limits

output_dir = "Results";

if save_files == true
    filename = sprintf('%s\\Relevant_curvature_profile.png',output_dir);
    saveas(fig, filename);
end

name_fig = sprintf('Relevant curvature profile: Lateral deviation');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[s]'); ylabel('[m]');
plot(tout,e1, 'LineWidth', 1)
axis normal

if save_files == true
    filename = sprintf('%s\\Relevant_curvature_profile_lateral_deviation.png',output_dir);
    saveas(fig, filename);
end


name_fig = sprintf('Relevant curvature profile: Vehicle trajectory');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[m]'); ylabel('[m]');
plot(Var_trajectory(:,1),Var_trajectory(:,2),'g',X,Y,'b--','LineWidth', 1.5);
axis normal
legend('Real','Ideal','Location', 'best')


if save_files == true
    filename = sprintf('%s\\Relevant_curvature_profile_trajectory.png',output_dir);
    saveas(fig, filename);
end

%% Skid-pad test

% TO DO

%% Obstacle avoidance manoeuvre

% TO DO

%%


%% Results plots

% TO DO

