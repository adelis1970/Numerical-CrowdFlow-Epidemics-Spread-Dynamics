
function y = FSG(phi, c, Dx, Dy)
% FSGuess - Fast Sweeping update for anisotropic Eikonal equation
% phi : current phi values
% c   : local speed function
% Dx  : grid spacing in x
% Dy  : grid spacing in y
%
% y   : updated phi after one full 4-sweep iteration

[Ny, Nx] = size(phi);
phiNEW = phi;
phiOLD = phi;

% Helper function: anisotropic Eikonal update at a single point
updatePoint = @(a,b,cij,Dx,Dy) ...
    solveQuadratic(a,b,cij,Dx,Dy);

% 4 sweeps
sweeps = {...
    [2, Ny-1, 1], [2, Nx-1, 1];  % left->right, top->bottom
    [2, Ny-1, 1], [Nx-1, 2, -1]; % right->left, top->bottom
    [Ny-1, 2, -1], [Nx-1, 2, -1];% right->left, bottom->top
    [Ny-1, 2, -1], [2, Nx-1, 1]  % left->right, bottom->top
};

for s = 1:4
    i_range = sweeps{s,1}(1):sweeps{s,1}(3):sweeps{s,1}(2);
    j_range = sweeps{s,2}(1):sweeps{s,2}(3):sweeps{s,2}(2);
    for i = i_range
        for j = j_range
            if phiOLD(i,j) ~= 0
                % Neighbor values
                a = min(phiOLD(max(i-1,1),j), phiOLD(min(i+1,Ny),j));
                b = min(phiOLD(i,max(j-1,1)), phiOLD(i,min(j+1,Nx)));
                
                % Check for one-sided update (if difference too large)
                if abs(a-b) >= c(i,j)*sqrt(Dx^2 + Dy^2)
                    phiNEW(i,j) = min(a + c(i,j)*Dy, b + c(i,j)*Dx);
                else
                    phiNEW(i,j) = updatePoint(a,b,c(i,j),Dx,Dy);
                end
                
                % Ensure monotonicity
                phiNEW(i,j) = min(phiNEW(i,j), phiOLD(i,j));
            end
        end
    end
    phiOLD = phiNEW; % update phiOLD after each sweep
end

y = phiNEW;

end

%-----------------------------
function phi_ij = solveQuadratic(a,b,cij,Dx,Dy)
% Solve the anisotropic quadratic:
% (phi - a)^2/Dy^2 + (phi - b)^2/Dx^2 = cij^2

% Coefficients
A = 1/Dy^2 + 1/Dx^2;
B = -2*(a/Dy^2 + b/Dx^2);
C = a^2/Dy^2 + b^2/Dx^2 - cij^2;

phi_ij = (-B + sqrt(B^2 - 4*A*C)) / (2*A);
end