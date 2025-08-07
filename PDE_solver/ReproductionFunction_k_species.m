function B = ReproductionFunction_k_species(a, lag, lag_sum, b, alpha, rhok) 
%% ReproductionFunction_k_species --------------------------------------
% purpose : Compute total births for a single species via
%           integral from tau to inf. rho(a) · n_i(a,t‑tau_i) · max{1‑alpha_i* sum_j n_j, 0} da.
%
% inputs  : a        – age grid (1×n)
%           n_lag    – lagged density of focal species (n×1)
%           lag_sum  – sum of lagged densities of all species (n×1)
%           b        – birth‑rate coefficient (scalar)
%           alpha    – competition coefficient (scalar)
%           rho      – patch‑age density (1×n)
%
% outputs : B        – scalar birth term for the focal species
%-----------------------------------------------------------------------




B = b*trapz(a,rhok'.*lag.*max(1-alpha*lag_sum,0)); % density dependent reproduction

end