clear all; 
close all;
clc

save_files = false;
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
%poles = [-1 -2 -10 -50]; % No integral contribution
poles = [-5 -7 -10 -15 -30];

% % No integral contribution
% B1_c = [0
%     CF/m
%     0
%     (CF*a)/Jz];

AM = cell(1,N_SAMPLES);
%K_lookup = zeros(N_SAMPLES,4); % No integral contribution
K_lookup = zeros(N_SAMPLES,5);
Kff_lookup = zeros(N_SAMPLES,1);

% No integral contribution
% for i = 1:N_SAMPLES 
%     V(i) = (DELTA_SPEED*i)/3.6; % speeds from 10 to 130 Km/h expressed in m/s
% 
%     AM{i} = [0 1 0 0
%     0 -(CF+CR)/(m*V(i)) (CF+CR)/m (CR*b-CF*a)/(m*V(i))
%     0 0 0 1
%     0 (CR*b-CF*a)/(Jz*V(i)) (CF*a-CR*b)/Jz -(CR*b^2+CF*a^2)/(Jz*V(i))];
% 
%     K_lookup(i,:) = place(AM{i},B1_c, poles);
%     Kff_lookup(i) = ((m*V(i)^2)/l)*(b/CF-a/CR+(a*K_lookup(i,3)/CR))+l-b*K_lookup(i,3);
% 
% end

B1_c = [0
    CF/m
    0
    (CF*a)/Jz
    0];

% Linear Quadratic Regulator
Q=diag([1 1 1 1 1]);
R=1;
P_lqr = zeros(N_SAMPLES,5);

for i = 1:N_SAMPLES 
    V(i) = (DELTA_SPEED*i)/3.6; % speeds from 10 to 130 Km/h expressed in m/s

    AM{i} = [0 1 0 0 0
    0 -(CF+CR)/(m*V(i)) (CF+CR)/m (CR*b-CF*a)/(m*V(i)) 0
    0 0 0 1 0
    0 (CR*b-CF*a)/(Jz*V(i)) (CF*a-CR*b)/Jz -(CR*b^2+CF*a^2)/(Jz*V(i)) 0
    1 0 0 0 0];
    
    % % Pole Placement
    K_lookup(i,:) = place(AM{i},B1_c, poles);
    Kff_lookup(i) = ((m*V(i)^2)/l)*(b/CF-a/CR+(a*K_lookup(i,3)/CR))+l-b*K_lookup(i,3);
    
    % Linear Quadratic Regulator
    % [K_lookup(i,:),~,P_lqr(i,:)]=lqr(AM{i}, B1_c, Q, R);
    % Kff_lookup(i) = ((m*V(i)^2)/l)*(b/CF-a/CR+(a*K_lookup(i,3)/CR))+l-b*K_lookup(i,3);

end

%% Gain variations in function of vehicle speed

vout = linspace(10,130,N_SAMPLES);

name_fig1 = sprintf('Gain variations in function of vehicle speed');
fig1 = figure('Name',name_fig1);
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('V'); ylabel('K');
plot(vout,K_lookup,'LineWidth', 1.5);
axis normal
legend('k1','k2','k3','k4','k_i','Location', 'best')

%% Analysis of single-track model with no control

velocities=[5:1:25,25:5:130]/3.6;
poles = zeros(length(velocities),2);

for i = 1:length(velocities)
    V = velocities(i);
    
    % State space definition
    A=[(-CF-CR)/(m*V),(-CF*a+CR*b-m*V^2)/(m*V^2);
        (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*V)];
    B=[CF/(m*V) CR/(m*V);
        (CF*a/Jz) -(CR*b/Jz)];
    C = [1,0
        0,1
        (-CR-CF)/(m*V^2),(-CF*a+CR*b)/(m*V^3)
        -1, -a/V
        -1, b/V
        (-CR-CF)/(m),(-CF*a+CR*b)/(m*V)];
    D = [0 0;
        0 0;
        CF/(m*V^2) CR/(m*V^2)
        1 0
        0 1
        CF/m CR/m];
    
    G = ss(A,B,C,D);
    [Wn,Z,P]=damp(G);
    
    poles(i,:)=P;
end

% poles
P1=poles(:,1);
P2=poles(:,2);
P1_Real=real(P1);
P1_Im=imag(P1);
P2_Real=real(P2);
P2_Im=imag(P2);

fig=figure('Name','Poles');
figure(fig)
scatter(P1_Real,P1_Im,[],1:1:length(velocities))
hold on
grid on
scatter(P2_Real,P2_Im,[],1:1:length(velocities))
xlabel('Real'),ylabel('Im'),title('Poles'),
set(gca,'FontName','Times New Roman','FontSize',12)

output_dir = "Results";

if save_files == true
    filename = sprintf('%s\\Single_track_model_eigenvalues.png',output_dir);
    saveas(fig, filename);
end

controlled = false;
Tsim = 6;
speed_profile = 1;
curvature_profile = 1;


sim("model.slx");

name_fig = sprintf('Step steering response of single-track model');
fig = figure('Name',name_fig);
subplot(2,1,1); % 2 rows, 1 column, first subplot
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[m]'); ylabel('[yaw_rate]');
plot(tout,beta,'b','LineWidth', 1.5);
axis normal

subplot(2,1,2); % 2 rows, 1 column, second subplot
hold on, grid on
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('[m]'); ylabel('[beta]');
plot(tout,yaw_rate,'b','LineWidth', 1.5);
axis normal

if save_files == true
    filename = sprintf('%s\\Step_steering_response_single_track_model.png',output_dir);
    saveas(fig, filename);
end

controlled = true;

%% Evaluation along the relevant curvature profile
%close all;

curvature_profile = 1;
speed_profile = 1;

% TO DO: evaluate for different tunings

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

%close all;

curvature_profile = 2;
speed_profile = 4;

% TO DO

Tsim = 100;

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
V = 80/3.6;

% TO DO

Tsim = 10;

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
legend('Real','Ideal','Location', 'best')


if save_files == true
    filename = sprintf('%s\\Avoidance_test_trajectory.png',output_dir);
    saveas(fig, filename);
end

