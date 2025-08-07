function F = koren_flux(u)
%% koren_flux -----------------------------------------------------------
% purpose : Return left‑biased interface states using the third‑order
%           Koren limiter (TVD) for a 1‑D upwind scheme.
%
% inputs  : u – cell‑centred values (n×1)
% outputs : F – reconstructed left states at cell interfaces (n×1)
%-----------------------------------------------------------------------
    n   = numel(u);
    F   = u;                           % donor values (boundaries)

    idx  = 3 : n-2;                    % interior cells
    fwd  = u(idx+1) - u(idx);          % Right hand flux
    back = u(idx)   - u(idx-1);        % Left hand flux

    % Slope ratio r = (u_i - u_{i-1}) / (u_{i+1} - u_i)
    r    = back ./ (fwd + (10e-6));        % eps avoids 0/0 in constant regions

    % ------- Koren limiter: phi(r) = max(0, min(2r, (1+2r)/3, 2)) ----------
    phi  = max(0, min( [ 2*r; (1+2*r)/3; 2*ones(size(r)) ] ));

    % Reconstructed left value at i+1/2
    F(idx) = u(idx) + 0.5 .* phi .* fwd;

end