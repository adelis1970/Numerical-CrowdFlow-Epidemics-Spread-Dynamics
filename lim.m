function r = lim(dp, dm, ~)
% Slope limiter for MUSCL (Van Albada + kappa form, after Toro)
% Inputs:
%   dp = forward difference  (R(i+1)-R(i))
%   dm = backward difference (R(i)-R(i-1))
%   ~  = unused (kept for compatibility)
% Output:
%   r  = limiter value

    kap = 0;        % can tune if needed
    tol = 1e-12;    % numerical safety

    % Replace tiny values to avoid division by zero
    if abs(dm) < tol
        dm = tol * sign(dm + tol);
    end
    if abs(dp) < tol
        dp = tol * sign(dp + tol);
    end

    r1 = dp / dm;

    if r1 >= 0
        denor = 1 - kap + (1 + kap)*r1;
        phir  = 2 / denor;

        % Van Albada limiter
        phi   = (r1*(1+r1)) / (1 + r1^2);

        % Conservative choice
        r = min(phi, phir);
    else
        r = 0;  % first-order if slopes of opposite sign
    end
end
