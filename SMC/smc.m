%% Sliding Mode Control Simulation
clear all; close all; clc;

% Simulation parameters
t_final = 10;       % Simulation time
Ts = 0.001;         % Sampling time
t = 0:Ts:t_final;   % Time vector

m = 1;              % Mass (kg)
b = 0.5;            % Damping coefficient (N·s/m)
k = 0.1;            % Spring constant (N/m)

% Reference signal (desired trajectory)
xd = sin(t);        % Position reference
dxd = cos(t);       % Velocity reference
ddxd = -sin(t);     % Acceleration reference

% Initial conditions
x = [0; 0];         % [position; velocity]

% SMC parameters
lambda = 5;         % Sliding surface slope (tuning parameter)
eta = 1.5;          % Switching gain (robustness against disturbances)
Phi = 0.1;          % Boundary layer thickness (reduces chattering)

% Preallocate for logging
x1_log = zeros(size(t));
x2_log = zeros(size(t));
u_log = zeros(size(t));
s_log = zeros(size(t));

% Main simulation loop
for i = 1:length(t)
    % Current states
    x1 = x(1);   %Position
    x2 = x(2);  %Velocity
    
    % Tracking error
    e = x1 - xd(i);
    de = x2 - dxd(i);
    
    
    % Sliding surface
    s = de + lambda*e;
    
    % System Dynamics 
    f_x = (-b*x2 - k*x1)/m; % Known dynamics (friction + spring)
    g_x = 1/m;  % Control input gain
    
    % Disturbance (unknown to the controller)
    d = 0.5*sin(2*t(i));  % Time-varying disturbance
    
    % Control Law (SMC)
    u_eq = -f_x + ddxd(i) - lambda*de;  % Equivalent control
    u_sw = -eta*sat(s/Phi);             % Switching control
    u = (1/g_x)*(u_eq + u_sw);
    
    % System Update (Euler Integration)
    dx2 = f_x + g_x*u + d;  % Acceleration (dynamics + control + disturbance)
    x(2) = x(2) + dx2*Ts;   % Update velocity
    x(1) = x(1) + x(2)*Ts;  % Update position
    
    % Log data
    x1_log(i) = x(1);  % Store position
    x2_log(i) = x(2);  % Store velocity
    u_log(i) = u;      % Store control input
    s_log(i) = s;      % Store sliding surface
end

%% result plotting
figure;
subplot(3,1,1);
plot(t, xd, 'r--', t, x1_log, 'b-', 'LineWidth', 1.5);
legend('Reference', 'Actual');
ylabel('Position');
title('Sliding Mode Control Performance');

subplot(3,1,2);
plot(t, s_log, 'LineWidth', 1.5);
ylabel('Sliding Surface');
title('Sliding Surface Behavior');

subplot(3,1,3);
plot(t, u_log, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Input');
title('Control Effort');

figure;
plot(x1_log - xd, x2_log - dxd, 'LineWidth', 1.5);
xlabel('Position Error');
ylabel('Velocity Error');
title('Phase Portrait of Tracking Errors');
grid on;

% Saturation function for boundary layer (must be at the end)
function out = sat(in)
    if abs(in) <= 1
        out = in;
    else
        out = sign(in);
    end
end