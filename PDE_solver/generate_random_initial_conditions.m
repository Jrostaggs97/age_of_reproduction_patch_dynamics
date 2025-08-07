function initcon_struct = generate_random_initial_conditions(k, tau_vec, seed)
%% generate_random_initial_conditions ----------------------------------
% purpose : Create k anonymous functions a↦n_i(a,0) with Gaussian bumps
%           centred near twice the maturation delay.
%
% inputs  : k       – number of species
%           tau_vec – maturation delays (k×1)
%           seed    – RNG seed (optional)
%
% outputs : initcon_struct – k×1 cell array of function handles
%-----------------------------------------------------------------------

    if nargin == 3
        rng(seed);  % Set seed for reproducibility
    end

    
    initcon_struct = cell(k, 1); % Preallocate cell array

    % Distribution parameters
    A_mu = 0.3; A_sd = 0.1;

    for i = 1:k
        A_i = normrnd(A_mu, A_sd);
        mu_i = normrnd(2*tau_vec(i), .2);

        % Create function handle
        initcon_struct{i} = @(a) A_i * exp(-((a - mu_i).^2) /(tau_vec(i)/10));
    end
end


