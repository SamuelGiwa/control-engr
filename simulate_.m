%% Non-Isothermal CSTR - Open Loop Simulation

%CSTR_model.slx for simulink implementation
clear; clc; close all;

%%  PARAMETERS 
V = 1.0;            % m^3
rho = 1000;         % kg/m^3
Cp = 4.18;          % kJ/kg-K

k0 = 8.46e6;        % 1/s
Ea = 65;            % kJ/mol
R  = 0.008314;      % kJ/mol-K
dH = -60;           % kJ/mol

UA = 20;            % kJ/s-K

CAin = 4.0;         % kmol/m^3
Tin  = 350;         % K

q  = 0.1;           % m^3/s (constant)
Tc = 300;           % K (constant)

%%  STEADY-STATE (for reference) 
ss_fun = @(x)[
    (q/V)*(CAin - x(1)) - k0*exp(-Ea/(R*x(2)))*x(1)
    (q/V)*(Tin - x(2)) + (-dH/(rho*Cp))*k0*exp(-Ea/(R*x(2)))*x(1) ...
    - (UA/(V*rho*Cp))*(x(2)-Tc)
];

x_ss = fsolve(ss_fun,[1;350],optimoptions('fsolve','Display','off'));
CA_ss = x_ss(1);
T_ss  = x_ss(2);

fprintf('Steady state:\nCA = %.3f kmol/m^3\nT  = %.2f K\n\n',CA_ss,T_ss);

%% = SIMULATION =
tspan = [0 200];
x0 = [CAin; Tin];   % initial condition

[t,x] = ode15s(@(t,x) cstr_ode(x), tspan, x0);

%% = PLOTTING =
figure;
subplot(2,1,1)
plot(t,x(:,1),'b','LineWidth',1.5)
yline(CA_ss,'k--','Steady state')
ylabel('C_A (kmol/m^3)')
grid on

subplot(2,1,2)
plot(t,x(:,2),'r','LineWidth',1.5)
yline(T_ss,'k--','Steady state')
ylabel('T (K)')
xlabel('Time (s)')
grid on

sgtitle('Open-Loop Non-Isothermal CSTR')

%%  ODE FUNCTION =
function dx = cstr_ode(x)
    CA = x(1);
    T  = x(2);

    % Parameters
    V = 1.0; rho = 1000; Cp = 4.18;
    k0 = 8.46e6; Ea = 65; R = 0.008314; dH = -60;
    UA = 20; CAin = 4; Tin = 350;
    q = 0.1; Tc = 300;

    % Reaction rate
    k = k0*exp(-Ea/(R*T));
    r = k*CA;

    % CSTR equations
    dCA = (q/V)*(CAin - CA) - r;
    dT  = (q/V)*(Tin - T) + (-dH/(rho*Cp))*r ...
          - (UA/(V*rho*Cp))*(T - Tc);

    dx = [dCA; dT];
end
