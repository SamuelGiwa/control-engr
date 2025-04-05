%% Kinetic Analysis of Pyrolysis using Model-Fitting Methods
% Implementation of Coats-Redfern and Criado's master plot methods

clc;
clear;
close all;

%% Sample Data Generation (Replace with actual experimental data)
T = linspace(300, 1000, 100)';       % Temperature (K)
alpha = linspace(0, 0.99, 100)';     % Conversion degree (0-1)
beta = 10/60;                        % Heating rate (K/s)
dalpha_dt = exp(-0.002*(T-400).^2);  % Sample mass loss rate (1/s)

% Universal gas constant
R = 8.314; % J/(mol·K)

%% Coats-Redfern Method Analysis
fprintf('=== Coats-Redfern Method Analysis ===\n');

% Define reaction models (g(alpha) functions)
models = {
    {'F1', @(a) -log(1 - a)},                % First-order
    {'F2', @(a) 1./(1 - a) - 1},             % Second-order
    {'D1', @(a) a.^2},                       % 1-D Diffusion
    {'D2', @(a) (1-a).*log(1-a) + a},        % 2-D Diffusion (Valensi)
    {'D2-AJ', @(a) ((1+a).^(1/3) - 1).^2},   % 2-D Anti-Jander
    {'D3', @(a) (1 - (1-a).^(1/3)).^2},     % 3-D Diffusion (Jander)
    {'D4', @(a) 1 - (2/3)*a - (1-a).^(2/3)} % 3-D Diffusion (G-B)
};

best_R2 = 0;
best_model = '';
best_Ea = 0;
best_A = 0;

for i = 1:length(models)
    model = models{i};
    name = model{1};
    g_func = model{2};
    
    try
        % Calculate left side of CR equation
        y = log(g_func(alpha)./(T.^2));
        
        % Prepare for linear regression (1/T vs y)
        x = 1./T;
        
        % Perform linear fit (y = mx + c)
        p = polyfit(x, y, 1);
        m = p(1);
        c = p(2);
        
        % Calculate R-squared
        y_pred = m*x + c;
        ss_res = sum((y - y_pred).^2);
        ss_tot = sum((y - mean(y)).^2);
        R2 = 1 - (ss_res / ss_tot);
        
        % Calculate kinetic parameters
        Ea = -m * R / 1000; % kJ/mol
        A = exp(c) * beta * Ea * 1000 / R;
        
        % Display results for current model
        fprintf('Model: %-6s | Ea: %6.1f kJ/mol | A: %8.2e 1/s | R2: %.4f\n', ...
                name, Ea, A, R2);
        
        % Check if this model fits better
        if R2 > best_R2
            best_R2 = R2;
            best_model = name;
            best_Ea = Ea;
            best_A = A;
        end
    catch
        fprintf('Error processing model %s\n', name);
    end
end

fprintf('\nBest fitting model: %s\n', best_model);
fprintf('Activation energy: %.1f kJ/mol\n', best_Ea);
fprintf('Pre-exponential factor: %.2e 1/s\n', best_A);
fprintf('R-squared: %.4f\n\n', best_R2);

%% Criado's Master Plot Method
fprintf('=== Criado''s Master Plot Method ===\n');

% Find index where alpha is closest to 0.5
[~, idx_05] = min(abs(alpha - 0.5));
T_05 = T(idx_05);
dadt_05 = dalpha_dt(idx_05);

% Calculate Z(alpha)/Z(0.5) for experimental data
Z_exp = (T./T_05).^2 .* (dalpha_dt./dadt_05);

% Theoretical models for Criado's plot
criado_models = {
    {'F1', @(a) 1-a, @(a) -log(1-a)}, % First-order
    {'F2', @(a) (1-a).^2, @(a) 1./(1-a) - 1}, % Second-order
    {'D1', @(a) 1./(2*a), @(a) a.^2}, % 1-D Diffusion
    {'D2', @(a) -1./log(1-a), @(a) (1-a).*log(1-a) + a}, % 2-D Valensi
    {'D2-AJ', @(a) 3*(1+a).^(2/3)./(2*((1+a).^(1/3)-1)), ... % 2-D Anti-Jander
              @(a) ((1+a).^(1/3)-1).^2},
    {'D3', @(a) 3*(1-a).^(2/3)./(2*(1-(1-a).^(1/3))), ... % 3-D Jander
              @(a) (1-(1-a).^(1/3)).^2}
};

% Generate theoretical curves
figure;
hold on;
grid on;

% Plot experimental data
plot(alpha, Z_exp, 'k-', 'LineWidth', 2, 'DisplayName', 'Experimental');

% Plot theoretical models
colors = lines(length(criado_models));
for i = 1:length(criado_models)
    model = criado_models{i};
    name = model{1};
    f = model{2};
    g = model{3};
    
    % Calculate Z(alpha)/Z(0.5)
    fg_alpha = f(alpha) .* g(alpha);
    fg_05 = f(0.5) * g(0.5);
    Z_ratio = fg_alpha ./ fg_05;
    
    plot(alpha, Z_ratio, '--', 'Color', colors(i,:), 'DisplayName', name);
end

xlabel('Conversion (\alpha)');
ylabel('Z(\alpha)/Z(0.5)');
title('Criado''s Master Plot');
legend('Location', 'best');
xlim([0 1]);
hold off;

%% Comparison of Results
fprintf('\n=== Comparison of Results ===\n');
fprintf('Both methods should ideally identify the same best-fitting model.\n');
fprintf('In this analysis:\n');
fprintf('Coats-Redfern best model: %s\n', best_model);
fprintf('Visual inspection of Criado''s plot should confirm this.\n');