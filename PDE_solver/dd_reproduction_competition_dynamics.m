function [time_grid, n_mat] = dd_reproduction_competition_dynamics(k, b_vec, mu_vec, tau_vec, alpha_vec, gamma, initcon_struct, Na, amax, tmax)
%-----------------------------------------------------------------------
% purpose : Solve a k‑species age‑structured PDE with maturation delay and
%           competition on reproductive output using dde23 in time and a flux‑
%           limited upwind scheme in age.
%
% inputs  : k            – number of species
%           b_vec        – birth‑rate coefficients           (k×1)
%           mu_vec       – mortality coefficients            (k×1)
%           tau_vec      – maturation delays (age of first reproduction)
%           alpha_vec    – LV competition coefficients       (k×1)
%           gamma        – patch‑age decay rate (ρ(a)=γe^{-γa})
%           initcon_struct –  k×1 cell of fns, one per species: a↦n(a,0)
%           Na           – number of age cells
%           amax         – maximum age in grid
%           tmax         – final simulation time
%
% outputs : time_grid    – 1×m vector of simulation times (from dde23)
%           n_mat        – (Na+1)·k × m matrix stacking species blocks
%
% notes   : – Spatial discretisation: donor‑cell + optional Koren/minmod
%           – Time integration: MATLAB's dde23 (explicit RK + interpolation)
%           – Assumes initial conditions passed in as a struct
%           – P in solver is total population. Each species age
%           distribution is stacked to work in solver environment
%-----------------------------------------------------------------------



%% Grid -----------------------------------------------------------------
da = amax/Na;                % uniform age step delta_a
a = 0:da:amax;                % age cells 

%% Pre‑compute delay indices -------------------------------------------
if isempty(tau_vec) == 1                % zero delay case
    tau_vec_1 = [0];                
    tau_a_idx_vec = int32(tau_vec_1./da)+1;
else
    tau_a_idx_vec = int32(tau_vec./da);   %idx where adults start
end 


%% Build initial condition matrix --------------------------------------
initcon_vec = zeros(length(a),k);
for i = 1:k
    initcon_vec(:,i) = initcon_struct{i,1}(a)';  % user define initial condition struct 
    initcon_vec(1:tau_a_idx_vec(i),i) = zeros(tau_a_idx_vec(i),1);  % enforce no adult population below age of maturity (a consistency condition)
end


%% Equilibrium patch‑age density ---------------------------------------
rho = gamma.*exp(-gamma.*a); % assumes steady state of McKendrick Von Foerster PDE with constant disturbance rate gamma

%% Call DDE solver ------------------------------------------------------
Z = tau_vec;                % vector of delays
hist_segment = initcon_vec; % repeat initial condition for history values 

% P stands for total population. 
sol = dde23(@(t,P,Z)...
     dd_reproduction_rhs_fun(t,P,Z,b_vec,mu_vec,alpha_vec,a,da,tau_vec,tau_a_idx_vec,rho,k),Z,hist_segment,[0,tmax]); 


%% Package outputs ------------------------------------------------------
time_grid = sol.x; % Extract time points from the solution
n_mat  = sol.y; % relabel for convenience 


end