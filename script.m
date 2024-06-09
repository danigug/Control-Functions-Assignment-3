clear all; 
close all;
clc

%% Parameters
m = 1575; % [Kg] vehicle_mass
Jz = 2875; % [kg*m^2] yaw_mass_moment_inertia
a = 1.3; % [m] front_semi_wheelbase
b = 1.5; % [m] rear_semi_wheelbase
l = a+b;
CF = 2*60000; % [N/rad] front_axle_cornering_stiffness
CR = 2*60000; % [N/rad] rear_axle_cornering_stiffness

delta_max = 25; % [deg] maximum_steering_angle

%% Look-up table of control gains
N_SAMPLES = 13;
MAX_SPEED = 130; % [km/h]
DELTA_SPEED = MAX_SPEED/N_SAMPLES;
poles = [-10 -20 -100 -200];

B1_c = [0
    CF/m
    0
    (CF*a)/Jz];

AM = cell(1,N_SAMPLES);
K_lookup = zeros(N_SAMPLES,4);
Kff_lookup = zeros(N_SAMPLES,1);

for i = 1:N_SAMPLES 
    V(i) = (DELTA_SPEED*i)/3.6; % speeds from 10 to 130 Km/h expressed in m/s

    AM{i} = [0 1 0 0
    0 -(CF+CR)/(m*V(i)) (CF+CR)/m (CR*b-CF*a)/(m*V(i))
    0 0 0 1
    0 (CR*b-CF*a)/(Jz*V(i)) (CF*a-CR*b)/Jz -(CR*b^2+CF*a^2)/(Jz*V(i))];

    K_lookup(i,:) = place(AM{i},B1_c, poles);
    Kff_lookup(i) = ((m*V(i)^2)/l)*(b/CF-a/CR+(a*K_lookup(i,3)/CR))+l-b*K_lookup(i,3);

end

%% Evaluation along the relevant curvature profile
close all;

curvature_profile = 1;
speed_profile = 1;

% TO DO: evaluate for different tunings

save_files = false;

V = 80/3.6;
Tsim = 200;

sim("model.slx");

name_fig = sprintf('Relevant curvature profile');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[s]'); ylabel('[m^-1]');
plot(tout,Kl_profile, 'LineWidth', 1)
axis tight;  
currentYlim = ylim;  
padding = 0.07 * (currentYlim(2) - currentYlim(1));  
ylim([currentYlim(1) - padding, currentYlim(2) + padding]);

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

close all;

curvature_profile = 2;
speed_profile = 4;

% TO DO

Tsim = 100;

save_files = true;
sim("model.slx");

name_fig = sprintf('Skidpad test: Lateral deviation');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[s]'); ylabel('[m]');
plot(tout,e1, 'LineWidth', 1)
axis normal

if save_files == true
    filename = sprintf('%s\\Skidpad_lateral_deviation.png',output_dir);
    saveas(fig, filename);
end

name_fig = sprintf('Skidpad: Reference and actual trajectory');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[m]'); ylabel('[m]');
plot(Var_trajectory(:,1),Var_trajectory(:,2),'g',X,Y,'b--','LineWidth', 1.5);
axis equal
legend('Real','Ideal','Location', 'best')


if save_files == true
    filename = sprintf('%s\\Skidpad_trajectory.png',output_dir);
    saveas(fig, filename);
end

%% Obstacle avoidance manoeuvre

close all;

curvature_profile = 4;
speed_profile = 1;
V = 85/3.6;

% TO DO

Tsim = 30;

save_files = true;
sim("model.slx");

name_fig = sprintf('Avoidance test: Lateral deviation');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[s]'); ylabel('[m]');
plot(tout,e1, 'LineWidth', 1)
axis normal

if save_files == true
    filename = sprintf('%s\\Avoidance_test_deviation.png',output_dir);
    saveas(fig, filename);
end

name_fig = sprintf('Avoidance test: Reference and actual trajectory');
fig = figure('Name',name_fig);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[m]'); ylabel('[m]');
plot(Var_trajectory(:,1),Var_trajectory(:,2),'g',X,Y,'b--','LineWidth', 1.5);
axis normal
legend('Real','Ideal','Location', 'best')


if save_files == true
    filename = sprintf('%s\\Avoidance_test_trajectory.png',output_dir);
    saveas(fig, filename);
end

%%


%% Results plots

% TO DO

