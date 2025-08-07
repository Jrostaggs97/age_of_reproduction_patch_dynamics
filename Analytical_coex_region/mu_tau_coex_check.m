function [S_i, S_j, mu_grid, tau_grid] = mu_tau_coex_check(alpha, gamma, b, ...
                                        mu_i_vec, tau_i_vec, ...
                                        mu_j, tau_j)
%-----------------------------------------------------------------------
% purpose : Compute invasion scores (S_i, S_j) on a grid where focal
%           species‑i varies in mortality μ_i and maturation age τ_i;
%           resident species‑j is fixed at (μ_j, τ_j).  Birth‑rate (b),
%           competition (α), and decay rate (γ) are shared.
%
% inputs  : alpha     – competition coefficient α            (scalar)
%           gamma     – decay rate γ                         (scalar)
%           b         – shared birth‑rate b                  (scalar)
%           mu_i_vec  – vector of focal mortalities μ_i      (rows)
%           tau_i_vec – vector of focal τ_i                  (cols)
%           mu_j      – resident mortality μ_j (scalar)
%           tau_j     – resident τ_j          (scalar)
%
% outputs : S_i       – |μ_i_vec| × |τ_i_vec| invasion scores for i
%           S_j       – corresponding matrix for resident j
%           mu_grid   – meshgrid of μ_i (R×C)
%           tau_grid  – meshgrid of τ_i (R×C)
%
% notes   : Handles four sign‑cases of P_i, P_j to mimic numerical
%           coexistence tests 
%-----------------------------------------------------------------------

%------------------ Grid construction ----------------------------------
mu_i_vec  = mu_i_vec(:);       % R rows
tau_i_vec = tau_i_vec(:)';     % C cols
[mu_i, tau_i] = meshgrid(mu_i_vec, tau_i_vec);  % C×R
mu_i  = mu_i';   % R×C
tau_i = tau_i';  % R×C

mu_grid = mu_i;
tau_grid = tau_i;

%------------------ Invader (spec i) terms -----------------------------
P_i = (b .* exp(-gamma .* tau_i)) ./ (gamma + mu_i) - 1;

A_i = (b.^2 .* gamma .* exp(-gamma .* tau_i)) ./ (mu_i.^2) .* ...
      (1/gamma - 2 ./ (gamma + mu_i) + 1 ./ (gamma + 2 .* mu_i));

%------------------ Resident (species‑j) terms --------------------------
P_j = (b * exp(-gamma .* tau_j)) ./ (gamma + mu_j) - 1;

A_j = (b^2 * gamma .* exp(-gamma .* tau_j)) ./ (mu_j.^2) .* ...
      (1/gamma - 2./(gamma + mu_j) + 1./(gamma + 2*mu_j));

delta = tau_i - tau_j;

B = gamma * (b^2) ./ (mu_i * mu_j) .* exp(-gamma * tau_i) .* ...
    (1/gamma ...
     - exp(-mu_j .* delta) ./ (gamma + mu_j) ...
     - 1 ./ (gamma + mu_i) ...
     + exp(-mu_j .* delta) ./ (gamma + mu_i + mu_j));


%------------------ Invasion scores ------------------------------------

% Preallocate S_i and S_j
S_i = zeros(size(P_i));
S_j = zeros(size(P_i));

% Logical masks
mask11 = (P_i > 0) & (P_j > 0);   % both positive
mask01 = (P_i <= 0) & (P_j > 0);  % only P_j positive
mask10 = (P_i > 0) & (P_j <= 0);  % only P_i positive
mask00 = (P_i <= 0) & (P_j <= 0); % both nonpositive

% Case: both P_i and P_j positive
S_i(mask11) = alpha .* A_j .* P_i(mask11) - alpha .* B(mask11) .* P_j;
S_j(mask11) = alpha .* A_i(mask11) .* P_j - alpha .* B(mask11) .* P_i(mask11);

% Case: only P_j positive
S_i(mask01) = -1;
S_j(mask01) =  1;

% Case: only P_i positive
S_i(mask10) =  1;
S_j(mask10) = -1;

% Case: both nonpositive
S_i(mask00) = -1;
S_j(mask00) = -1;




end
