function F = minmod_flux(u)

%% minmod_flux ----------------------------------------------------------
% purpose : Second‑order TVD minmod limiter reconstruction.
%
% inputs  : u – cell‑centred values (n×1)
% outputs : F – reconstructed left states at interfaces (n×1)
%-----------------------------------------------------------------------

n  = numel(u);          % total points
F  = u;                 % doner values

idx       = 3:n-2;       % interior indices that need the correction
forward   = u(idx+1) - u(idx);     % Right hand flux
backward  = u(idx)   - u(idx-1);   % Left hand flux

% Reconstructed left value at interface (i + 1/2) with left handed bias
% minmod limiter, phi(r) = max(0, min(1,r))
F(idx) = u(idx) + 0.5 .* ...
         0.5 .* (sign(forward)+sign(backward)) .* min(abs(forward), abs(backward));
end