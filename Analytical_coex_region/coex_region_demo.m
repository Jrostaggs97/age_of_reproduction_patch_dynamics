%-----------------------------------------------------------------------
% purpose : Example script that computes and plots coexistence outcomes
%           for a mu–tau trade-off using *mu_tau_coex_check.m*.
%
% inputs  : (edit here if you want different parameter sweeps)
%           alpha     – competition coefficient                (scalar)
%           gamma     – patch-age decay rate                   (scalar)
%           b         – per-capita birth rate (shared)         (scalar)
%           mu_j      – resident adult mortality               (scalar)
%           tau_j     – resident age of first reproduction     (scalar)
%           mu_i_vec  – vector of invader mortality values     (rows)
%           tau_i_vec – vector of invader τ values             (cols)
%
% outputs : - Figure: pcolor plot labelled "Coexistence region, μ–τ trade-off"
%           - coex   : R×C matrix of outcome codes
%                      0 = no winner, 1 = resident wins,
%                      2 = invader wins, 3 = coexistence
% Notes:    - typically scale parameters to be on 100 year scalings (e.g.
%             tau = .3 -> means age of first reproduction is 30 years and
%             mu = 7 means a per capita mortality rate of 7 adult individuals per
%             100 years)
%-----------------------------------------------------------------------

clc; clear;

%% Parameters -----------------------------------------------------------
alpha = 0.1;
gamma = 0.5;
b     = 24;

mu_j  = 7;        % resident mortality
tau_j = 0.3;      % resident age of first reproduction (tau)

delta_mu  = 0.005;
delta_tau = 0.005;

mu_i_vec  = 2:delta_mu: mu_j;   % rows in the plot
tau_i_vec = tau_j:delta_tau: 2.7;  % cols in the plot

%% Compute S-scores and classify outcomes ------------------------------
[S_i, S_j, mu_grid, tau_grid] = ...
    mu_tau_coex_check(alpha, gamma, b, ...
                      mu_i_vec, tau_i_vec, ...
                      mu_j,    tau_j);

i_check = 2*(S_i > 0);     % 2 if invader succeeds, 0 otherwise
j_check =     (S_j > 0);   % 1 if resident succeeds, 0 otherwise
coex    = i_check + j_check;

% Anchor color scale by setting four corner cells (never shown)
coex(1,1:4) = [0 1 2 3];

%% Plot -----------------------------------------------------------------
figure;
p = pcolor(tau_grid, mu_grid, coex);
set(p, 'EdgeColor', 'none');
ylim([min(mu_i_vec)+delta_mu, max(mu_i_vec)]);

% Custom 4-color map: white | red | blue | black
colormap([1 1 1; 1 0 0; 0 0 1; 0 0 0]);

% Colorbar with explicit labels
c = colorbar;
c.Ticks      = [.375 .9375 1.6875 2.625];
c.TickLabels = {'No winner', 'Resident wins', ...
                'Invader wins', 'Coexistence'};
c.Label.String = 'Outcome';
c.FontSize     = 10;

xlabel('Invader age of first reproduction, $\tau_i$', ...
       'FontSize', 25, 'Interpreter', 'latex');
ylabel('Invader adult mortality, $\mu_i$', ...
       'FontSize', 25, 'Interpreter', 'latex');
title('Coexistence region, $\mu$--$\tau$ trade-off', ...
      'FontSize', 25, 'Interpreter', 'latex');
