function F = Approx_Roe_PW(Rr,Rl,Ur,Ul,Vr,Vl,M)
% Approximate Riemann solver (Roe by default; optional Rusanov).
% M==0 -> x-face; M==1 -> y-face

    c0    = 1.2;        % "sound speed" in your ped model
    eps11 = 1.0e-4;     % small number used in your code
    tiny  = 1.0e-12;    % numerical guard

    % Roe averages with guards to avoid NaNs
    Rr_pos = max(Rr,0);  Rl_pos = max(Rl,0);
    sr = sqrt(Rr_pos);   sl = sqrt(Rl_pos);
    den = sr + sl;

    if den > eps11
        u = (Ur*sr + Ul*sl) / den;
        v = (Vr*sr + Vl*sl) / den;
    else
        u = 0; v = 0;
    end

    % eigenvalues for x and y systems (diagonal matrices written explicitly)
    eigx = [c0+u  0     0;
            0     u     0;
            0     0  -c0+u];

    eigy = [c0+v  0     0;
            0     v     0;
            0     0  -c0+v];

    % Choose orientation and assemble flux pieces
    if M==0
        % right/left fluxes in x-direction
        F_R = [Rr*Ur;             Rr*Ur^2 + c0^2*Rr;  Rr*Ur*Vr];
        F_L = [Rl*Ul;             Rl*Ul^2 + c0^2*Rl;  Rl*Ul*Vl];
        q_r = [Rr;                Rr*Ur;              Rr*Vr];
        q_l = [Rl;                Rl*Ul;              Rl*Vl];

        % right-eigenvector matrix you used (keep as-is)
        RR  = [1     0     1;
               c0+u  0   -c0+u;
               v     1     v];

        % Harten entropy fix on the three eigenvalues
        lam  = [eigx(1,1); eigx(2,2); eigx(3,3)];
        epsv = [max([0, lam(1)-(c0+Ul), (c0+Ur)-lam(1)]);
                max([0, lam(2)-(Ul)   , (Ur)    -lam(2)]);
                max([0, lam(3)-(-c0+Ul),(-c0+Ur)-lam(3)])];
        Lam = zeros(3,1);
        for k=1:3
            if abs(lam(k)) >= epsv(k)
                Lam(k) = abs(lam(k));
            else
                ek = max(epsv(k), tiny);  % guard
                Lam(k) = (lam(k)^2 + ek^2)/(2*ek);
            end
        end
        eig_f = diag(Lam);

    else
        % right/left fluxes in y-direction
        F_R = [Rr*Vr;             Rr*Ur*Vr;           Rr*Vr^2 + c0^2*Rr];
        F_L = [Rl*Vl;             Rl*Ul*Vl;           Rl*Vl^2 + c0^2*Rl];
        q_r = [Rr;                Rr*Ur;              Rr*Vr];
        q_l = [Rl;                Rl*Ul;              Rl*Vl];

        % right-eigenvector matrix you used for y (keep as-is)
        RR  = [1     0     1;
               u     1     u;
               c0+v  0   -c0+v];

        % Harten entropy fix  on the three eigenvalues
        lam  = [eigy(1,1); eigy(2,2); eigy(3,3)];
        epsv = [max([0, lam(1)-(c0+Vl), (c0+Vr)-lam(1)]);
                max([0, lam(2)-(Vl)   , (Vr)    -lam(2)]);
                max([0, lam(3)-(-c0+Vl),(-c0+Vr)-lam(3)])];
        Lam = zeros(3,1);
        for k=1:3
            if abs(lam(k)) >= epsv(k)
                Lam(k) = abs(lam(k));
            else
                ek = max(epsv(k), tiny);
                Lam(k) = (lam(k)^2 + ek^2)/(2*ek);
            end
        end
        eig_f = diag(Lam);
    end

    % Roe  matrix |A| ≈ R |Λ| R^{-1}
    
    Qa = (RR * abs(eig_f)) / RR;

    %Flux: Roe
    F_roe = 0.5*(F_R + F_L) - 0.5 * Qa * (q_r - q_l);

    %Optional: Rusanov (LLF) switch (set to true to use)
    use_rusanov = false;  % <— flip to true if you want Rusanov flux for testing
    if use_rusanov
        if M==0
            smax = max([abs(Ur)+c0, abs(Ul)+c0]);   % x-direction speeds
        else
            smax = max([abs(Vr)+c0, abs(Vl)+c0]);   % y-direction speeds
        end
        F_Rus = 0.5*(F_R + F_L) - 0.5*smax*(q_r - q_l);
        F = F_Rus;
    else
        F = F_roe;
    end
end
