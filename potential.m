
% Computes a steady velocity airflow field UG=[ug, vg] for the
% drif-reaction-diffusion equation for the infection potenrial β

%To be executed first if vendilation is included in the simulation scenario

clear; 

% Domain setup (as in the main program)
Lx = 10; Ly = 10;
Nx = 101; Ny = 101;
dx = Lx/(Nx-1); dy = Ly/(Ny-1);

x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

%%Initialize potential & source
psi = zeros(Ny, Nx);
src = zeros(Ny, Nx); % Define source term if needed (not used in this version)

% Poisson solver parameters
maxIter = 5000;

% Custom boundary regions besed on current discritization
N11=41; N12=61; %for 2m ventilation ducts 

uin = 10; %vendilation speed at inflow
uout = -uin;
%% Iterative solver (finite difference)
for iter = 1:maxIter
    % Vectorized interior update
    psi_new = psi;
    psi_new(2:end-1,2:end-1) = ( ...
        (psi(3:end,2:end-1) + psi(1:end-2,2:end-1))*dx^2 + ...
        (psi(2:end-1,3:end) + psi(2:end-1,1:end-2))*dy^2 + ...
        dx^2 * dy^2 .* src(2:end-1,2:end-1) ) / (2*(dx^2 + dy^2));
    
    % Neumann boundary conditions
    psi_new(1,:) = psi_new(2,:);        % bottom
    psi_new(end,:) = psi_new(end-1,:);  % top
    psi_new(:,1) = psi_new(:,2);        % left
    psi_new(:,end) = psi_new(:,end-1);  % right
    
    % Custom flux for selected left/right boundary region for ventilation inflow and
    % outblow ducts
    
    psi_new(N11:N12,1) = psi_new(N11:N12,2) + dx*uin;
    psi_new(N11:N12,end) = psi_new(N11:N12,end-1) - dx*uin;
    
    % Update psi
    psi = psi_new;
end

% Compute airflow velocity field
[ug, vg] = gradient(-psi, dx, dy);

% Plot potential
figure;
contourf(x, y, psi, 50, 'LineStyle', 'none');
colorbar;
xlabel('x'); ylabel('y');
title('Potential Field');

% Plot velocity field and streamlines
figure;
quiver(x(1:2:end), y(1:2:end), ug(1:2:end,1:2:end), vg(1:2:end,1:2:end), 4);
xlabel('x'); ylabel('y');
hold on;
h = streamslice(x, y, ug, vg);
set(h,'color','red','LineWidth',1);
axis([0 Lx 0 Ly]);
hold off;