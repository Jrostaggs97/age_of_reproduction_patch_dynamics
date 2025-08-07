function [S_i,S_j,b_grid,tau_grid] = b_tau_coex_check( ...
                       alpha,gamma,mu, ...          % scalars
                       b_j,tau_j, ...               % resident scalars
                       b_i_vec,tau_i_vec)           % vectors (rows=b_i, cols=tau_i)
%-----------------------------------------------------------------------
% purpose : Compute invasion scores (S_i, S_j) on a grid where the focal
%           species‑i varies in birth‑rate (b_i) and maturation age (tau_i),
%           while the resident species‑j is fixed at (b_j, tau_j).  Both
%           species share common mortality (mu), competition (alpha), and
%           disturbance decay (gamma).
%
% inputs  : alpha   – Lotka–Volterra competition coefficient (scalar)
%           gamma   – patch‑age decay rate gamma (scalar)
%           mu      – mortality rate mu   (scalar)
%           b_j     – resident birth‑rate               (scalar)
%           tau_j   – resident age at first reproduction (scalar)
%           b_i_vec – vector of focal birth‑rates  (rows in output)
%           tau_i_vec – vector of focal tau_i        (cols in output)
%
% outputs : S_i      – |b_i_vec| × |τ_i_vec| matrix of invasion scores for i
%           S_j      – corresponding matrix for resident j
%           b_grid   – meshgrid of b_i (R×C)
%           tau_grid – meshgrid of tau_i (R×C)
%
% notes   : Positive S_k ⇒ species‑k can invade when rare.  See analytical
%           derivation in Appendix of the manuscript.
%-----------------------------------------------------------------------


%------------------ Grid construction ----------------------------------
% force orientations
b_i_vec   = b_i_vec(:);        % R rows  (birth-rate axis)
tau_i_vec = tau_i_vec(:)';     % C cols  (age axis)
[b_grid,tau_grid] = meshgrid(b_i_vec,tau_i_vec);  % C×R
b_grid   = b_grid';   % R×C
tau_grid = tau_grid'; % R×C

%------------------ Resident (species‑j) terms --------------------------
P_j = (b_j * exp(-gamma * tau_j)) / (gamma + mu) - 1;
A_j = (b_j^2 * gamma * exp(-gamma * tau_j)) / mu^2 * ...
      (1/gamma - 2/(gamma + mu) + 1/(gamma + 2*mu));

%------------------ Invader (species‑i) grid terms ------------------------
P_i = (b_grid .* exp(-gamma .* tau_grid)) ./ (gamma + mu) - 1;

A_i = (b_grid.^2 .* gamma .* exp(-gamma .* tau_grid)) ./ mu^2 .* ...
      (1/gamma - 2./(gamma + mu) + 1./(gamma + 2*mu));

delta = tau_grid - tau_j;

B = gamma * (b_grid .* b_j) ./ mu^2 .* exp(-gamma .* tau_grid) .* ...
    ( 1/gamma ...
    - exp(-mu .* delta) ./ (gamma + mu) ...
    - 1/(gamma + mu) ...
    + exp(-mu .* delta) ./ (gamma + 2*mu) );

%------------------ Invasion scores ------------------------------------
S_i =  alpha * A_j .* P_i - alpha * B .* P_j;
S_j =  alpha * A_i .* P_j - alpha * B .* P_i;
end
