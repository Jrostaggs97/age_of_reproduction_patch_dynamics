%-----------------------------------------------------------------------
% purpose : Example simulation of a two-species PDE with a μ–τ trade-off.
%           Generates random Gaussian initial conditions, integrates with
%           *dd_reproduction_competition_dynamics*, and plots total
%           abundance N_i(t) for each species.
%
% inputs  : (edit these blocks to explore different scenarios)
%           • Global parameters
%                 gamma       – patch-age decay rate                (scalar)
%                 tmax        – final simulation time               (scalar)
%           • Species vectors (k = 2 here)
%                 b_vec       – per-capita birth rates              (1×k)
%                 mu_vec      – adult mortality rates               (1×k)
%                 tau_vec     – ages of first reproduction          (1×k)
%                 alpha_vec   – LV competition coefficients         (1×k)
%           • Discretisation
%                 delat_a      – desired age-step (sets Na)
%
% outputs : • time_grid – 1×m vector of output times
%           • n_mat     – (k·Na+1)×m matrix of age distributions
%           • Figure    – Total abundance curves N_i(t)
%-----------------------------------------------------------------------

clc; clear;
seed = 4;
rng(seed);                                % reproducible random seed

%% Global parameters ----------------------------------------------------
gamma = 0.5;                           % patch-age decay
tmax  = 25;                           % simulation horizon

%% Species parameters ---------------------------------------------------
k        = 2;
b_vec    = [24,  24];
mu_vec   = [5.89, 7.00];
tau_vec  = [1.633, 0.300];
alpha_vec = [0.1, 0.1];

%% Spatial grid ---------------------------------------------------------
amax = max( round(-(1/gamma)*log(1e-8/gamma)), 20 );
da   = 0.01;
Na   = ceil(amax/da);                  % integer number of cells
a    = linspace(0, amax, Na+1);        % age grid (row vector)

%% Initial conditions ---------------------------------------------------
initcon_struct = generate_random_initial_conditions(k, tau_vec, 4);
                                           % cell of anonymous fns

%% Solve PDE system -----------------------------------------------------
[time_grid, n_mat] = dd_reproduction_competition_dynamics( ...
    k, b_vec, mu_vec, tau_vec, alpha_vec, ...
    gamma, initcon_struct, Na, amax, tmax);

%% Compute total abundance N_i(t) ---------------------------------------
rho   = gamma * exp(-gamma * a);       % patch-age density
n_len = numel(a);                      % length of one species block
N     = zeros(k, numel(time_grid));    % pre-allocate

for i = 1:k
    for t = 1:length(time_grid)
        N(i,t) = trapz(a, rho'.*n_mat((i-1)*length(a)+1:i*length(a),t));
    end
end


%% Total Abundance Plot -----------------------------------------------------------------
figure;
plot(time_grid, N', 'LineWidth', 1.5);   % transpose → lines = species
xlabel('time, $t$',           'FontSize', 25, 'Interpreter', 'latex');
ylabel('$N(t)$',            'FontSize', 25, 'Interpreter', 'latex');
title('Two-species total abundance, $\gamma = 0.5$', ...
      'FontSize', 25, 'Interpreter', 'latex');

legend(cellfun(@(x) sprintf('$\\tau = %.3f$', x), num2cell(tau_vec), ...
               'UniformOutput', false), ...
       'FontSize', 14, 'Interpreter', 'latex');


%% Final age-distribution snapshot -------------------------------------
figure;
for i = 1:k
    rows = (i-1)*n_len + (1:n_len);     % rows for species i
    plot(a, n_mat(rows, end), 'LineWidth', 1.5);   % final-time column
    hold on
end
xlabel('patch age, $a$', 'FontSize', 25, 'Interpreter', 'latex');
ylabel('$n(a,t_{\mathrm{final}})$', ...
       'FontSize', 25, 'Interpreter', 'latex');
title('Final age distributions', ...
      'FontSize', 25, 'Interpreter', 'latex');

legend(cellfun(@(x) sprintf('$\\tau = %.3f$', x), num2cell(tau_vec), ...
               'UniformOutput', false), ...
       'FontSize', 14, 'Interpreter', 'latex');

