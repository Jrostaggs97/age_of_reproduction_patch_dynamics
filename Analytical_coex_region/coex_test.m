clc;
clear;
alpha = .1;
gamma = .5;
b= 40; %

mu_j = 7;
tau_j = .3;

delta_tau = .005;
%delta_b = .005;
delta_mu = .005;

% a quicker/cleaner way to select for given tau and mu would be nice
tau_i_vec = [tau_j:delta_tau:2.7];
%b_i_vec = [2:delta_b:mu_j];%
mu_i_vec = [2:delta_mu:mu_j]; %... 
%mu_i_vec = [5:delta_mu:25];
[S_i, S_j, mu_grid, tau_grid] = mu_tau_coex_check(alpha, gamma, b,mu_i_vec, tau_i_vec, mu_j, tau_j);

i_check = 2*(S_i >0);
j_check = S_j >0;
coex = i_check + j_check; 

coex(1,1) = 0; %use as anchor point for color We need to anchor every color
coex(1,2) = 1; 
coex(1,3) = 2;
coex(1,4) = 3;
% 
% [S_i_int, S_j_int] = integral_coex_check(gamma, tau_i_vec(1), b, mu_i_vec(1), alpha, tau_j, b, mu_j, alpha);
% 
% S_i - S_i_int
% S_j - S_j_int
% --- Plot ---
figure;
% Use pcolor for plotting
p = pcolor(tau_grid, mu_grid, coex);
set(p, 'EdgeColor', 'none');
ylim([min(mu_i_vec)+delta_mu,max(mu_i_vec)])


% Custom colormap: white, red, blue, black
cmap = [1 1 1;      % 0 = no winner
        1 0 0;      % 1 = resident wins
        0 0 1;      % 2 = invader wins
        0 0 0];     % 3 = coexistence
colormap(cmap);

% Colorbar with custom labels
c = colorbar;
c.Ticks = [.75/2, mean([.75,1.5]), mean([1.5,2.25]),  mean([2.25,3])];
c.TickLabels = {'No winner', 'Resident wins', 'Invader wins', 'Coexistence'};
c.Label.String = 'Outcome';
c.FontSize = 10;


xlabel('Invader age of first reproduction, $\tau_i$',"Fontsize", 25, 'Interpreter', 'latex');
ylabel('Invader adult mortality, $\mu_i$', "Fontsize", 25,'Interpreter', 'latex');
title("Coexistence region, $\mu-\tau$ trade off","Fontsize", 25, 'Interpreter', 'latex');
%% b tau trade off 
% alpha = .1;
% gamma = .5;
% mu = 7;
% 
% b_j = 12;
% tau_j = .3;%
% 
% 
% 
% delta_tau = .005;
% delta_b = .005;
% 
% % a quicker/cleaner way to select for given tau and mu would be nice
% tau_i_vec = [tau_j:delta_tau:2.7];
% b_i_vec = [b_j:delta_b:28];%[.01:delta_mu:mu_j]; %... 
% %mu_i_vec = [5:delta_mu:25];
% [S_i, S_j, b_grid, tau_grid] = b_tau_coex_check(alpha,gamma,mu, ...
%                                      b_j,tau_j, ...
%                                      b_i_vec,tau_i_vec);
% 
% i_check = 2*(S_i >0);
% j_check = S_j >0;
% coex = i_check + j_check; 
% 
% coex(1,1) = 0; %use as anchor point for color We need to anchor every color
% coex(1,2) = 1; 
% coex(1,3) = 2;
% coex(1,4) = 3;
% % 
% % [S_i_int, S_j_int] = integral_coex_check(gamma, tau_i_vec(1), b, mu_i_vec(1), alpha, tau_j, b, mu_j, alpha);
% % 
% % S_i - S_i_int
% % S_j - S_j_int
% % --- Plot ---
% figure;
% % Use pcolor for plotting
% p = pcolor(tau_grid, b_grid, coex);
% set(p, 'EdgeColor', 'none');
% ylim([min(b_i_vec)+delta_b,max(b_i_vec)])
% 
% 
% % Custom colormap: white, red, blue, black
% cmap = [1 1 1;      % 0 = no winner
%         1 0 0;      % 1 = resident wins
%         0 0 1;      % 2 = invader wins
%         0 0 0];     % 3 = coexistence
% colormap(cmap);
% 
% % Colorbar with custom labels
% c = colorbar;
% c.Ticks = [.75/2, mean([.75,1.5]), mean([1.5,2.25]),  mean([2.25,3])];
% % c.TickLabels = {'No winner', 'Resident wins', 'Invader wins', 'Coexistence'};
% % c.Label.String = 'Outcome';
% c.FontSize = 10;
% 
% 
% xlabel('Invader age of first reproduction, $\tau_i$',"Fontsize", 25, 'Interpreter', 'latex');
% ylabel('Invader recruitment to adulthood, $b_i$', "Fontsize", 25,'Interpreter', 'latex');
% title("Coexistence region, $b-\tau$ trade off","Fontsize", 25, 'Interpreter', 'latex');
%%  mu-tau random sampling

%title("Coexistence region, $\mu-\tau$ trade off $(\alpha, \gamma, b) = (" +num2str(alpha) + ", " +num2str(gamma)+", " + num2str(b)+")$","Fontsize", 25, 'Interpreter', 'latex');
% 
% % Exclude anchor point row from sampling
% coex_sub = coex(2:end, :);
% mu_grid_sub = mu_grid(2:end, :);
% tau_grid_sub = tau_grid(2:end, :);
% 
% % Find linear indices of each class
% idx_resident_wins = find(coex_sub == 1);
% idx_invader_wins  = find(coex_sub == 2);
% idx_coex          = find(coex_sub == 3);
% 
% 
% % Helper to randomly sample 10 indices, avoiding empty cases
% sample_safe = @(idx) idx(randsample(numel(idx), min(10, numel(idx))));
% 
% s_resident = sample_safe(idx_resident_wins);
% s_invader  = sample_safe(idx_invader_wins);
% s_coex     = sample_safe(idx_coex);
% 
% % Combine all samples
% all_samples = [s_resident; s_invader; s_coex];
% labels = [repmat("j wins", length(s_resident), 1); ...
%           repmat("i wins", length(s_invader), 1); ...
%           repmat("coex",   length(s_coex), 1)];
% 
% % Convert linear indices to subscripts
% [row_idx, col_idx] = ind2sub(size(coex_sub), all_samples);
% 
% % Construct output struct
% for k = 1:length(all_samples)
%     mu_i = mu_grid_sub(row_idx(k), col_idx(k));
%     tau_i = tau_grid_sub(row_idx(k), col_idx(k));
% 
%     samples(k).gamma = gamma;
%     samples(k).alpha_vec = [alpha, alpha];
%     samples(k).b_vec = [b, b];
%     samples(k).mu_vec = [mu_i, mu_j];
%     samples(k).tau_vec = [tau_i, tau_j];
%     samples(k).outcome = labels(k);
% end
% 
% %Optional: save to file
% save('coex_samples_1.mat', 'samples');
% 
% % Save coexistence matrix
% coex_matrix = coex; % rename for clarity
% save('coex_matrix_1.mat', 'coex_matrix', 'mu_grid', 'tau_grid');

%% mu tau grid sampling
% Define the grid resolution
n_pts = 10;  % You can change this to 50, 100, etc.

% Exclude anchor point (row 1) from sampling
coex_sub     = coex(2:end, :);
mu_grid_sub  = mu_grid(2:end, :);
tau_grid_sub = tau_grid(2:end, :);

% Get coordinate range
mu_vals  = mu_grid_sub(:);
tau_vals = tau_grid_sub(:);

mu_min  = min(mu_vals);
mu_max  = max(mu_vals);
tau_min = min(tau_vals);
tau_max = max(tau_vals);

% Generate uniform grid in mu and tau space
[mu_i_grid, tau_i_grid] = meshgrid(linspace(mu_min, mu_max, n_pts), ...
                                   linspace(tau_min, tau_max, n_pts));

% Convert to vectors for looping
b_i_vec  = mu_i_grid(:);
tau_i_vec = tau_i_grid(:);

% Interpolate to nearest coex matrix values (to assign outcome)
% We'll use griddedInterpolant-style logic to find closest grid index
mu_axis  = mu_grid(2:end, 1);     % assuming rows vary by mu
tau_axis = tau_grid(2, :);        % assuming cols vary by tau

% Preallocate sample struct array
samples(n_pts^2, 1) = struct();

for k = 1:numel(b_i_vec)
    % Get current test pair
    mu_i  = b_i_vec(k);
    tau_i = tau_i_vec(k);

    % Find closest grid points in original coex matrix
    [~, mu_idx]  = min(abs(mu_axis - mu_i));
    [~, tau_idx] = min(abs(tau_axis - tau_i));

    % Fetch outcome from coex matrix
    coex_val = coex(mu_idx + 1, tau_idx);  % +1 to re-align with original coex

    % Assign label
    switch coex_val
        case 1
            label = "j wins";
        case 2
            label = "i wins";
        case 3
            label = "coex";
        otherwise
            label = "unknown";
    end

    % Build sample entry
    samples(k).gamma      = gamma;
    samples(k).alpha_vec  = [alpha, alpha];
    samples(k).b_vec      = [b, b];
    samples(k).mu_vec     = [mu_i, mu_j];
    samples(k).tau_vec    = [tau_i, tau_j];
    samples(k).outcome    = label;
end


hold on;

% === Step 2: Overlay green circles at sampled (tau_i, mu_i) ===
tau_samples = arrayfun(@(s) s.tau_vec(1), samples);
mu_samples  = arrayfun(@(s) s.mu_vec(1),  samples);

plot(tau_samples, mu_samples, 'go', 'MarkerSize', 6, 'LineWidth', 1.5);

% Save output
save('coex_samples_grid.mat', 'samples');

% Also save grid and matrix for visualization
coex_matrix = coex;save('coex_matrix_1.mat', 'coex_matrix', 'mu_grid', 'tau_grid');
