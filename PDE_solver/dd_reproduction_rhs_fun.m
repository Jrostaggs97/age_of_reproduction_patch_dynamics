function [v] = dd_reproduction_rhs_fun(t,P,Z,b_vec,mu_vec,alpha_vec,a,da,tau_vec,tau_a_idx_vec,rho,k)
%% dd_reproduction_rhs_fun ----------------------------------------------
% purpose : Right‑hand‑side callback for *dde23*; returns dP/dt for the
%           stacked population vector.
%
% inputs  : t      – current time (scalar, unused but required by dde23)
%           P      – current state vector ((Na+1)·k ×1)
%           Z      – matrix of delayed states, Z(:,i) = state at t‑τ_i
%           b_vec  – birth rates (k×1) | mu_vec – deaths (k×1)
%           alpha_vec – LV competition coefficients (k×1)
%           a, da  – age grid and step | tau_vec – delay vector
%           tau_a_idx_vec – integer indices where adulthood starts
%           rho    – patch‑age density | k – number of species
%
% outputs : dPdt   – time derivative, same size as P

% notes   : – Age flux can be either Koren or minmod limiter
%           – Minmod is second order [in smooth regions] and symmetric
%           – Koren is third order [in smooth regions] NOT symmetric
%           (symmetry is not so important in this case since wave velocity
%           (age) is strictly positive. (namely characteristic of 1) )
%-----------------------------------------------------------------------

    n = numel(a);              % number of age classes
    v_mat = zeros(n, k);       % each column i will be dP_i/dt over ages

    for i = 1:k
        idx    = (i-1)*n + (1:n);      % slice for species i
        spec_i = P(idx);               % current pop of species i
        blockZ = Z(:, i);               % lagged pop of species i
        lag_all = reshape(blockZ, [n, k]);     % reshape once to sum across spp  
        lag_i   = lag_all(:, i);      % species i’s lagged pop
        
        spec_sum_i = sum(lag_all, 2);      %delayed sum for species interactions: sum_j n_j(a,t‑tau_j)
        
        % ---------- Demographic terms -------------------------------------
        Repo_i  = ReproductionFunction_k_species(a, lag_i, spec_sum_i, b_vec(i), alpha_vec(i), rho); 
       
        %Flux_i  = minmod_flux(spec_i);
        Flux_i  = koren_flux(spec_i); % numerical flux (can swap to minmod)

        % Finite‑volume update ---------------------------------------------
        tau_idx = tau_a_idx_vec(i);  % assemble age of maturation indices
        n_i     = zeros(n,1);        % pre allocate vector of zeros 

        n_i(2:tau_idx) = -(1/da) * (Flux_i(2:tau_idx) - Flux_i(1:tau_idx-1)); %update for juveniles (just flux)

        n_i(tau_idx+1:n) = Repo_i - mu_vec(i).* spec_i(tau_idx+1:n) - (1/da) * (Flux_i(tau_idx+1:n) - Flux_i(tau_idx:n-1)); %update for adults (uses all demographic processes)

        
        v_mat(:,i) = n_i;          % store result
    end

    v = reshape(v_mat, [], 1);       %Flatten back into a (k*n)×1 column for dde23

end