%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second-order model solved by Roe's scheme for pedestrians' flow
% and Rusanov for the others models using MUSCL reconstruction
%for crowd flow towards an exit 2m wide
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; %if vendilation is included run first potential.m and comment this

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions & parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_t = 0.01;                   %Time step as to satisfy the CFL condition
D_x = 0.1; D_y = 0.1;         %Discretization parameters dx and dy
tow = 0.6; vf = 1.4; Rm = 6;  % Relaxation time, free flow speed and max density
maxx = 10; maxy = 10;         %Romm space dimensions

exit_length = 2; % Exit length placed at the center of right boundary

HV = 1e12;      % High value for the FSM
S = 1.2e-3;     % Diffusion coefficient for β

lose = 0.5;     % for particle losses in the drift--diffusion-reaction equation
kap = 0;        % for MUSCL reconstruction

i0 = 0.04;      % infection rate
eps1 = 1e-4; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 0:D_x:maxx;
y = 0:D_y:maxy;
[X, Y] = meshgrid(x, y);

Nx = length(x);
Ny = length(y);

% Number of cells needed
n_cells_exit = round(exit_length / D_y);

% Center index
center_idx = round(Ny / 2);

% Compute start and end indices (inclusive) for exit
NEXY1 = center_idx - floor((n_cells_exit-1)/2);
NEXY2 = center_idx + ceil((n_cells_exit-1)/2);

fprintf('Exit indices: %d to %d\n', NEXY1, NEXY2);

% Initialize fields
R_new = zeros(Ny, Nx);
U_new = zeros(Ny, Nx);
V_new = zeros(Ny, Nx);

% Initialize diffusion/other fields
b0 = zeros(Ny, Nx);
b1 = b0;
b = b0;

SUM_EXP_EXIT1 = 0;
SUM_EXP_EXIT2 = 0;
SUM_EXP_EXIT = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions for density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask_R = (X >= 1 & X <= 5) & (Y >= 2.5 & Y <= 7.5);
R_new(mask_R) = 2.4;

% Initial velocity fields
U_new(:) = 0;
V_new(:) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Infected, susceptible and masked percentage of pedestrians 
% (can be changed as needed for different scenarios)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
irow_ic = 0.25 * R_new;   % infected percentage
srow_ic = 0.75 * R_new;   % susceptible percentage
vrow_ic = 0.0 * R_new;   % vaccinated/masked pedestrians
erow = zeros(size(R_new)); % exposed initially zero

% Copy for simulation updates
irow = irow_ic;
srow = srow_ic;
vrow = vrow_ic;

iirow = irow;
ssrow = srow;
vvrow = vrow;
eerow = erow;

% Velocity copies
u = U_new;
v = V_new;

% Previous step storage
Rin = R_new;
Uin = U_new;
Vin = V_new;


% Create a movie for evolved quantities and a results file 
%
writerObj = VideoWriter('Test_Case', 'mpeg-4'); %to create an evolution movie
writerObj.FrameRate = 50;
open(writerObj);

fileID = fopen('RESULTS.m','w'); % to record time pedestrians and exposed densities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:7200 %time loop
fprintf('t_spep = %d time = %f\n',t,t*D_t); 

clf;
%figure(1)
subplot(3,1,1);
surfc(x,y, R_new); colorbar; xlim([0,maxx]); ylim([0,maxy]); view (0,90) %;clim([0,4]);
shading INTERP
 title('Total Pedestrians Density $\rho$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
 xlabel('$x(m)$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
 ylabel('$y(m)$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
 set(gca,'FontName','Times New Roman','FontSize',20)
 
 subplot(3,1,2);
 surfc(x,y,b(1:Ny,1:Nx));colorbar; xlim([0,maxx]); ylim([0,maxy]); view (0,90) %;clim([0,0.03]);
 shading INTERP
 title('Infection Coefficient $\beta_I$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
 xlabel('$x(m)$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
 ylabel('$y(m)$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
 set(gca,'FontName','Times New Roman','FontSize',20)

subplot(3,1,3);
surfc(x,y,eerow);colorbar; xlim([0,maxx]); ylim([0,maxy]); view (0,90)%;clim([0,0.3]);
%quiver(x,y,U_new,V_new); xlim([0,50]); ylim([0,20]);
colormap(jet)
shading INTERP
title('Exposed Pedestrians Density $\rho^E$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
xlabel('$x(m)$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
ylabel('$y(m)$','interpreter', 'latex','FontSize',18,'FontName','Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',20)

sgtitle(sprintf('Time: %.1f s', t*D_t),'FontName','Times New Roman','FontSize',20);

 pause (0.001);
 Frame = getframe(gcf) ;
 writeVideo(writerObj, Frame);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Indices for interior cells
i=2:Nx+1;
j=2:Ny+1;

% Copy new solution into main grid
R(j,i)=R_new(1:Ny,1:Nx);
U(j,i)=U_new(1:Ny,1:Nx);
V(j,i)=V_new(1:Ny,1:Nx);

% -----------------------------
% Ghost Cell Boundary Conditions
% -----------------------------

% R (density) - zero flux at walls
R(1,i)    = R(2,i);       % bottom
R(Ny+2,i) = R(Ny+1,i);    % top
R(j,1)    = R(j,2);       % left
R(j,Nx+2) = R(j,Nx+1);    % right (except exit, see below)

U(1,i)    = U(2,i);  
U(Ny+2,i) = U(Ny+1,i); 
U(j,1)    = -U(j,2); 
U(j,Nx+2) = -U(j,Nx+1);

V(1,i)    =  -V(2,i);
V(Ny+2,i) =  -V(Ny+1,i);
V(j,1)    =  V(j,2);
V(j,Nx+2) =  V(j,Nx+1);

% Set exit boundary

R(NEXY1:NEXY2,Nx+2:Nx+2) = R(NEXY1:NEXY2,Nx+1:Nx+1);
U(NEXY1:NEXY2,Nx+2:Nx+2)  = vf;
V(NEXY1:NEXY2,Nx+2:Nx+2)  = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 UTMP = U; %Store pedestrians' velocities to be used for the SEIS model
 VTMP = V;
 RTMP= R;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize and solve the eikonal equation 

% Cost function (1/speed field)

speed = UV(R);              % desired speed from density
speed(speed < 1e-6) = 1e-6; % safeguard against division by zero
cost  = 1 ./ speed;         % eikonal "cost" field

% Initialize potential field

phi = HV * ones(Ny+2, Nx+2);

%Right exit, set phi=0 at boundary

phi(NEXY1:NEXY2, Nx+2) = 0;

%  Fast Sweeping Method iterations

tol    = 1e-9;
change = inf;

while change > tol
    phi_old = phi;
    phi     = FSG(phi, cost, D_x, D_y);  % Zhao (2005) sweeps
    change  = max(abs(phi(:) - phi_old(:)));
end

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Growd Flow Module 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

clear R_new U_new V_new Qold f_s FL FR GL GR ...
g_s Qold_04 Qold_34 Q_new Q_new_04 Q_new_34 ...
R0 U0 V0 R_34 U_34 V_34

i_int = 2:Nx+1;
j_int = 2:Ny+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Crowd Flow Module - One Time Step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y-SWEEP

for j = j_int
    for i = i_int
       
        % MUSCL reconstruction
        
        if i==2 || i==Nx+1 || j==2 || j==Ny+1 %first-order at the boundary
            FL = Approx_Roe_PW(R(j,i),R(j,i-1),U(j,i),U(j,i-1),V(j,i),V(j,i-1),0);
            FR = Approx_Roe_PW(R(j,i+1),R(j,i),U(j,i+1),U(j,i),V(j,i+1),V(j,i),0);
        else
            % Left reconstruction
            RL = R(j,i-1) + 0.25*lim(R(j,i)-R(j,i-1),R(j,i-1)-R(j,i-2),1)*((1-kap)*(R(j,i-1)-R(j,i-2)) + (1+kap)*(R(j,i)-R(j,i-1)));
            UL = U(j,i-1) + 0.25*lim(U(j,i)-U(j,i-1),U(j,i-1)-U(j,i-2),1)*((1-kap)*(U(j,i-1)-U(j,i-2)) + (1+kap)*(U(j,i)-U(j,i-1)));
            VL = V(j,i-1) + 0.25*lim(V(j,i)-V(j,i-1),V(j,i-1)-V(j,i-2),1)*((1-kap)*(V(j,i-1)-V(j,i-2)) + (1+kap)*(V(j,i)-V(j,i-1)));

            RR = R(j,i) - 0.25*lim(R(j,i+1)-R(j,i),R(j,i)-R(j,i-1),1)*((1-kap)*(R(j,i+1)-R(j,i)) + (1+kap)*(R(j,i)-R(j,i-1)));
            UR = U(j,i) - 0.25*lim(U(j,i+1)-U(j,i),U(j,i)-U(j,i-1),1)*((1-kap)*(U(j,i+1)-U(j,i)) + (1+kap)*(U(j,i)-U(j,i-1)));
            VR = V(j,i) - 0.25*lim(V(j,i+1)-V(j,i),V(j,i)-V(j,i-1),1)*((1-kap)*(V(j,i+1)-V(j,i)) + (1+kap)*(V(j,i)-V(j,i-1)));

            FL = Approx_Roe_PW(RR,RL,UR,UL,VR,VL,0);

            % Right reconstruction
            RL = R(j,i) + 0.25*lim(R(j,i+1)-R(j,i),R(j,i)-R(j,i-1),1)*((1-kap)*(R(j,i)-R(j,i-1)) + (1+kap)*(R(j,i+1)-R(j,i)));
            UL = U(j,i) + 0.25*lim(U(j,i+1)-U(j,i),U(j,i)-U(j,i-1),1)*((1-kap)*(U(j,i)-U(j,i-1)) + (1+kap)*(U(j,i+1)-U(j,i)));
            VL = V(j,i) + 0.25*lim(V(j,i+1)-V(j,i),V(j,i)-V(j,i-1),1)*((1-kap)*(V(j,i)-V(j,i-1)) + (1+kap)*(V(j,i+1)-V(j,i)));

            RR = R(j,i+1) - 0.25*lim(R(j,i+2)-R(j,i+1),R(j,i+1)-R(j,i),1)*((1-kap)*(R(j,i+2)-R(j,i+1)) + (1+kap)*(R(j,i+1)-R(j,i)));
            UR = U(j,i+1) - 0.25*lim(U(j,i+2)-U(j,i+1),U(j,i+1)-U(j,i),1)*((1-kap)*(U(j,i+2)-U(j,i+1)) + (1+kap)*(U(j,i+1)-U(j,i)));
            VR = V(j,i+1) - 0.25*lim(V(j,i+2)-V(j,i+1),V(j,i+1)-V(j,i),1)*((1-kap)*(V(j,i+2)-V(j,i+1)) + (1+kap)*(V(j,i+1)-V(j,i)));

            FR = Approx_Roe_PW(RR,RL,UR,UL,VR,VL,0);
        end

        % Update solution along x
       
        Qold = [R(j,i); R(j,i)*U(j,i); R(j,i)*V(j,i)];
        q_new = Qold - D_t/D_x * (FR - FL);

        % Source term along x

        r = q_new(1);
        u = 2*q_new(2)*q_new(1)/(q_new(1)^2 + max(q_new(1)^2,1e-8));
        UR = vf*exp(-7.5*(r/Rm)^2)*sign(u);
        if i==Nx+1
            UR2 = UR^2*(phi(j,i)-phi(j,i-1))/D_x;
        elseif i==2
            UR2 = UR^2*(phi(j,i+1)-phi(j,i))/D_x;
        else
            UR2 = UR^2*(phi(j,i+1)-phi(j,i-1))/(2*D_x);
        end
        f_s = [0; (-UR2*r - q_new(2))/tow; 0];
        Q_new = q_new + D_t*f_s;

        R0(j,i) = Q_new(1);
        U0(j,i) = 2*Q_new(2)*Q_new(1)/(Q_new(1)^2 + max(Q_new(1)^2,1e-8)); %Normalize velocities
        V0(j,i) = 2*Q_new(3)*Q_new(1)/(Q_new(1)^2 + max(Q_new(1)^2,1e-8));

        % Enforce small values threshold
        if R0(j,i) <= eps
            R0(j,i) = 0; U0(j,i) = 0; V0(j,i) = 0;
        end
        if abs(U0(j,i)) <= eps, U0(j,i) = 0; end
        if abs(V0(j,i)) <= eps, V0(j,i) = 0; end
    end
end


% Boundary update after x-sweep

i = 2:Nx+1; j = 2:Ny+1;
R0(1,i) = R0(2,i); R0(Ny+2,i) = R0(Ny+1,i);
R0(j,1) = R0(j,2); R0(j,Nx+2) = R0(j,Nx+1);
U0(1,i) = U0(2,i); U0(Ny+2,i) = U0(Ny+1,i);
U0(j,1) = -U0(j,2); U0(j,Nx+2) = -U0(j,Nx+1);
V0(1,i) = -V0(2,i); V0(Ny+2,i) = -V0(Ny+1,i);
V0(j,1) = V0(j,2); V0(j,Nx+2) = V0(j,Nx+1);

% Outflow boundary (right exit)
R0(NEXY1:NEXY2,Nx+2) = R0(NEXY1:NEXY2,Nx+1);
U0(NEXY1:NEXY2,Nx+2) = vf;
V0(NEXY1:NEXY2,Nx+2) = 0;


% y-SWEEP 

for j = j_int
    for i = i_int
       
        % MUSCL reconstruction
       
        if i==2 || i==Nx+1 || j==2 || j==Ny+1
            GL = Approx_Roe_PW(R0(j,i),R0(j-1,i),U0(j,i),U0(j-1,i),V0(j,i),V0(j-1,i),1);
            GR = Approx_Roe_PW(R0(j+1,i),R0(j,i),U0(j+1,i),U0(j,i),V0(j+1,i),V0(j,i),1);
        else
            % Bottom reconstruction
            RL = R0(j-1,i) + 0.25*lim(R0(j,i)-R0(j-1,i),R0(j-1,i)-R0(j-2,i),1)*((1-kap)*(R0(j-1,i)-R0(j-2,i)) + (1+kap)*(R0(j,i)-R0(j-1,i)));
            UL = U0(j-1,i) + 0.25*lim(U0(j,i)-U0(j-1,i),U0(j-1,i)-U0(j-2,i),1)*((1-kap)*(U0(j-1,i)-U0(j-2,i)) + (1+kap)*(U0(j,i)-U0(j-1,i)));
            VL = V0(j-1,i) + 0.25*lim(V0(j,i)-V0(j-1,i),V0(j-1,i)-V0(j-2,i),1)*((1-kap)*(V0(j-1,i)-V0(j-2,i)) + (1+kap)*(V0(j,i)-V0(j-1,i)));

            RR = R0(j,i) - 0.25*lim(R0(j+1,i)-R0(j,i),R0(j,i)-R0(j-1,i),1)*((1-kap)*(R0(j+1,i)-R0(j,i)) + (1+kap)*(R0(j,i)-R0(j-1,i)));
            UR = U0(j,i) - 0.25*lim(U0(j+1,i)-U0(j,i),U0(j,i)-U0(j-1,i),1)*((1-kap)*(U0(j+1,i)-U0(j,i)) + (1+kap)*(U0(j,i)-U0(j-1,i)));
            VR = V0(j,i) - 0.25*lim(V0(j+1,i)-V0(j,i),V0(j,i)-V0(j-1,i),1)*((1-kap)*(V0(j+1,i)-V0(j,i)) + (1+kap)*(V0(j,i)-V0(j-1,i)));

            GL = Approx_Roe_PW(RR,RL,UR,UL,VR,VL,1);

            % Top reconstruction
            RL = R0(j,i) + 0.25*lim(R0(j+1,i)-R0(j,i),R0(j,i)-R0(j-1,i),1)*((1-kap)*(R0(j,i)-R0(j-1,i)) + (1+kap)*(R0(j+1,i)-R0(j,i)));
            UL = U0(j,i) + 0.25*lim(U0(j+1,i)-U0(j,i),U0(j,i)-U0(j-1,i),1)*((1-kap)*(U0(j,i)-U0(j-1,i)) + (1+kap)*(U0(j+1,i)-U0(j,i)));
            VL = V0(j,i) + 0.25*lim(V0(j+1,i)-V0(j,i),V0(j,i)-V0(j-1,i),1)*((1-kap)*(V0(j,i)-V0(j-1,i)) + (1+kap)*(V0(j+1,i)-V0(j,i)));

            RR = R0(j+1,i) - 0.25*lim(R0(j+2,i)-R0(j+1,i),R0(j+1,i)-R0(j,i),1)*((1-kap)*(R0(j+2,i)-R0(j+1,i)) + (1+kap)*(R0(j+1,i)-R0(j,i)));
            UR = U0(j+1,i) - 0.25*lim(U0(j+2,i)-U0(j+1,i),U0(j+1,i)-U0(j,i),1)*((1-kap)*(U0(j+2,i)-U0(j+1,i)) + (1+kap)*(U0(j+1,i)-U0(j,i)));
            VR = V0(j+1,i) - 0.25*lim(V0(j+2,i)-V0(j+1,i),V0(j+1,i)-V0(j,i),1)*((1-kap)*(V0(j+2,i)-V0(j+1,i)) + (1+kap)*(V0(j+1,i)-V0(j,i)));

            GR = Approx_Roe_PW(RR,RL,UR,UL,VR,VL,1);
        end

        % Update solution along y

        Qold = [R0(j,i); R0(j,i)*U0(j,i); R0(j,i)*V0(j,i)];
        q_new = Qold - D_t/D_y * (GR - GL);

        % Source term along y

        r = q_new(1);
        v = 2*q_new(3)*q_new(1)/(q_new(1)^2 + max(q_new(1)^2,1e-8));
        VR = vf*exp(-7.5*(r/Rm)^2)*sign(v);

        if j==Ny+1
            VR2 = VR^2*(phi(j,i)-phi(j-1,i))/D_y;
        elseif j==2
            VR2 = VR^2*(phi(j+1,i)-phi(j,i))/D_y;
        else
            VR2 = VR^2*(phi(j+1,i)-phi(j-1,i))/(2*D_y);
        end

        g_s = [0;0;(-VR2*r - q_new(3))/tow];
        Q_new = q_new + D_t*g_s;

        R_new(j-1,i-1) = Q_new(1);
        U_new(j-1,i-1) = 2*Q_new(2)*Q_new(1)/(Q_new(1)^2 + max(Q_new(1)^2,1e-8));
        V_new(j-1,i-1) = 2*Q_new(3)*Q_new(1)/(Q_new(1)^2 + max(Q_new(1)^2,1e-8));

        % Threshold
        if R_new(j-1,i-1) <= eps
            R_new(j-1,i-1)=0; U_new(j-1,i-1)=0; V_new(j-1,i-1)=0;
        end
        if abs(U_new(j-1,i-1)) <= eps, U_new(j-1,i-1)=0; end
        if abs(V_new(j-1,i-1)) <= eps, V_new(j-1,i-1)=0; end
    end
end


 ug=zeros(Ny,Nx); %Air flow velocities set to if no ventilation is imposed
 vg=zeros(Ny,Nx); %else comments this and run potential.m first

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Infection potential based on drift-diffusionreaction equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Rusanov numerical fluxes
    F_x = zeros(Ny, Nx); 
    F_y = zeros(Ny, Nx);

    % Compute numerical fluxes in the x-direction

    for i = 1:Nx-1
        for j = 1:Ny

          if any(i<=2) || any(i>=Nx-2)

         max_abs_u_x = max(abs(ug(j,i)), abs(ug(j,i+1)));
         F_x(j,i) = 0.5*(b0(j,i)*ug(j,i)+ b0(j,i+1)*ug(j,i+1))- 0.5*max_abs_u_x*(b0(j,i+1) - b0(j,i));

         else
       ugr = ug(j,i+1) - 0.25*lim(ug(j,i+2)-ug(j,i+1),ug(j,i+1)-ug(j,i),1)*((1-kap)*(ug(j,i+2)-ug(j,i+1))+(1+kap)*(ug(j,i+1)-ug(j,i)));
       ugl = ug(j,i)   + 0.25*lim(ug(j,i+1)-ug(j,i),ug(j,i)-ug(j,i-1),1)*((1-kap)*(ug(j,i)-ug(j,i-1))+(1+kap)*(ug(j,i+1)-ug(j,i)));

       b0r = b0(j,i+1) - 0.25*lim(b0(j,i+2)-b0(j,i+1),b0(j,i+1)-b0(j,i),1)*((1-kap)*(b0(j,i+2)-b0(j,i+1))+(1+kap)*(b0(j,i+1)-b0(j,i)));
       b0l = b0(j,i)   + 0.25*lim(b0(j,i+1)-b0(j,i),b0(j,i)-b0(j,i-1),1)*((1-kap)*(b0(j,i)-b0(j,i-1))+(1+kap)*(b0(j,i+1)-b0(j,i)));

       max_abs_u_x = max(abs(ugr), abs(ugl));
        F_x(j,i) = 0.5*(b0r*ugr+ b0l*ugl)- 0.5*max_abs_u_x*(b0r - b0l);
            end
        end
    end

 % Compute numerical fluxes in the y-direction
    for i = 1:Nx
        for j = 1:Ny-1

            if any(j<=2) || any(j>=Ny-2)
            max_abs_v_y = max(abs(vg(j,i)), abs(vg(j+1,i)));
            F_y(j,i) = 0.5*(b0(j,i)*vg(j,i)+ b0(j+1,i)*vg(j+1,i))- 0.5*max_abs_v_y*(b0(j+1,i) - b0(j,i));  

            else
           vgr = vg(j+1,i) - 0.25*lim(vg(j+2,i)-vg(j+1,i),vg(j+1,i)-vg(j,i),1)*((1-kap)*(vg(j+2,i)-vg(j+1,i))+(1+kap)*(vg(j+1,i)-vg(j,i)));
           vgl = vg(j,i)   + 0.25*lim(vg(j+1,i)-vg(j,i),vg(j,i)-vg(j-1,i),1)*((1-kap)*(vg(j,i)-vg(j-1,i))+(1+kap)*(vg(j+1,i)-vg(j,i)));
           
           b0r = b0(j+1,i) - 0.25*lim(b0(j+2,i)-b0(j+1,i),b0(j+1,i)-b0(j,i),1)*((1-kap)*(b0(j+2,i)-b0(j+1,i))+(1+kap)*(b0(j+1,i)-b0(j,i)));
           b0l = b0(j,i)   + 0.25*lim(b0(j+1,i)-b0(j,i),b0(j,i)-b0(j-1,i),1)*((1-kap)*(b0(j,i)-b0(j-1,i))+(1+kap)*(b0(j+1,i)-b0(j,i)));

           max_abs_v_y = max(abs(vgr), abs(vgl));
           F_y(j,i) = 0.5*(b0r*vgr+ b0l*vgl)- 0.5*max_abs_v_y*(b0r - b0l);
            end
            
        end
    end

% Compute diffusion terms

    diff_x = S*(b0(2:end-1, 3:end) - 2 * b0(2:end-1, 2:end-1) + b0(2:end-1, 1:end-2)) / D_x^2;
    diff_y = S*(b0(3:end, 2:end-1) - 2 * b0(2:end-1, 2:end-1) + b0(1:end-2, 2:end-1)) / D_y^2;

    tmp_src =zeros(Ny,Nx);
% Time integration
    RT = iirow + ssrow + eerow + vvrow;
    tmp_src(2:end-1, 2:end-1)= iirow(2:end-1, 2:end-1)./(RT(2:end-1, 2:end-1)+eps);
  for j=2:Ny-1
      for i=2:Nx-1
          if RT(j,i) <= 10^(-3)
             tmp_src(j,i) =0;
        end
      end
  end

    b1(2:end-1, 2:end-1) = b0(2:end-1, 2:end-1) ...
        - D_t * (F_x(2:end-1, 2:end-1) - F_x(2:end-1, 1:end-2))/D_x ...
        - D_t * (F_y(2:end-1, 2:end-1) - F_y(1:end-2, 2:end-1))/D_y ...
        + D_t*(diff_x + diff_y-lose*b0(2:end-1, 2:end-1)+ tmp_src(2:end-1, 2:end-1));

  % Apply zero-flux boundary conditions
    b1(:, 1) = b1(:, 2); % Left boundary
    b1(:, end) =b1(:, end-1); % Right boundary
    b1(1, :) = b1(2, :); % Bottom boundary
    b1(end, :) =b1(end-1, :); % Top boundary

% Update for the next iteration
     b0 = b1;
     b=i0*b0;% the actual infection coefficient

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Solve SEIS model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 u=UTMP; %pedestrian velocities from the crowd flow model
 v=VTMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
i=2:Nx+1;
j=2:Ny+1;

irow(j,i) = iirow(1:Ny,1:Nx);
srow(j,i) = ssrow(1:Ny,1:Nx);
erow(j,i) = eerow(1:Ny,1:Nx);
vrow(j,i) = vvrow(1:Ny,1:Nx);


irow(1,i)  = irow(2,i); %Down
irow(Ny+2,i) = irow(Ny+1,i); %top
irow(j,1)  = irow(j,2);        %L
irow(j,Nx+2) = irow(j,Nx+1); 

srow(1,i)  = srow(2,i); %Down
srow(Ny+2,i) = srow(Ny+1,i); %top
srow(j,1)  = srow(j,2);        %L
srow(j,Nx+2) = srow(j,Nx+1);

erow(1,i)  = erow(2,i); %Down
erow(Ny+2,i) = erow(Ny+1,i); %top
erow(j,1)  = erow(j,2);        %L
erow(j,Nx+2) = erow(j,Nx+1);

vrow(1,i)  = vrow(2,i); %Down
vrow(Ny+2,i) = vrow(Ny+1,i); %top
vrow(j,1)  = vrow(j,2);        %L
vrow(j,Nx+2) = vrow(j,Nx+1);


% Compute Rusanov numerical fluxes
    Fi_x = zeros(Ny+2, Nx+1); % Initialize flux array
    Fi_y = zeros(Ny+1, Nx+2); 
    Fs_x = zeros(Ny+2, Nx+1); 
    Fs_y = zeros(Ny+1, Nx+2); 
    Fe_x = zeros(Ny+2, Nx+1); 
    Fe_y = zeros(Ny+1, Nx+2);
    Fv_x = zeros(Ny+2, Nx+1); 
    Fv_y = zeros(Ny+1, Nx+2); 

    % Compute numerical fluxes in the x-direction using MUSCL recons.

    for j = 1:Ny+2
        for i = 1:Nx+1
      
      if any(i<=2) || any(i >=Nx)

      max_abs_u_x = max(abs(u(j,i)), abs(u(j,i+1)));
      Fi_x(j,i) = 0.5*(irow(j,i)*u(j,i)+irow(j,i+1)*u(j,i+1))- 0.5*max_abs_u_x*(irow(j,i+1)-irow(j,i));
      Fs_x(j,i) = 0.5*(srow(j,i)*u(j,i)+srow(j,i+1)*u(j,i+1))- 0.5*max_abs_u_x*(srow(j,i+1)-srow(j,i));
      Fe_x(j,i) = 0.5*(erow(j,i)*u(j,i)+erow(j,i+1)*u(j,i+1))- 0.5*max_abs_u_x*(erow(j,i+1)-erow(j,i));
      Fv_x(j,i) = 0.5*(vrow(j,i)*u(j,i)+vrow(j,i+1)*u(j,i+1))- 0.5*max_abs_u_x*(vrow(j,i+1)-vrow(j,i));
      else
      ur = u(j,i+1) - 0.25*lim(u(j,i+2)-u(j,i+1),u(j,i+1)-u(j,i),1)*((1-kap)*(u(j,i+2)-u(j,i+1))+(1+kap)*(u(j,i+1)-u(j,i)));
      ul = u(j,i)   + 0.25*lim(u(j,i+1)-u(j,i),u(j,i)-u(j,i-1),1)*((1-kap)*(u(j,i)-u(j,i-1))+(1+kap)*(u(j,i+1)-u(j,i)));

      irowr = irow(j,i+1) - 0.25*lim(irow(j,i+2)-irow(j,i+1),irow(j,i+1)-irow(j,i),1)*((1-kap)*(irow(j,i+2)-irow(j,i+1))+(1+kap)*(irow(j,i+1)-irow(j,i)));
      irowl = irow(j,i)   + 0.25*lim(irow(j,i+1)-irow(j,i),irow(j,i)-irow(j,i-1),1)*((1-kap)*(irow(j,i)-irow(j,i-1))+(1+kap)*(irow(j,i+1)-irow(j,i)));

      srowr = srow(j,i+1) - 0.25*lim(srow(j,i+2)-srow(j,i+1),srow(j,i+1)-srow(j,i),1)*((1-kap)*(srow(j,i+2)-srow(j,i+1))+(1+kap)*(srow(j,i+1)-srow(j,i)));
      srowl = srow(j,i)   + 0.25*lim(srow(j,i+1)-srow(j,i),srow(j,i)-srow(j,i-1),1)*((1-kap)*(srow(j,i)-srow(j,i-1))+(1+kap)*(srow(j,i+1)-srow(j,i)));
   
      erowr = erow(j,i+1) - 0.25*lim(erow(j,i+2)-erow(j,i+1),erow(j,i+1)-erow(j,i),1)*((1-kap)*(erow(j,i+2)-erow(j,i+1))+(1+kap)*(erow(j,i+1)-erow(j,i)));
      erowl = erow(j,i)   + 0.25*lim(erow(j,i+1)-erow(j,i),erow(j,i)-erow(j,i-1),1)*((1-kap)*(erow(j,i)-erow(j,i-1))+(1+kap)*(erow(j,i+1)-erow(j,i)));

      vrowr = vrow(j,i+1) - 0.25*lim(vrow(j,i+2)-vrow(j,i+1),vrow(j,i+1)-vrow(j,i),1)*((1-kap)*(vrow(j,i+2)-vrow(j,i+1))+(1+kap)*(vrow(j,i+1)-vrow(j,i)));
      vrowl = vrow(j,i)   + 0.25*lim(vrow(j,i+1)-vrow(j,i),vrow(j,i)-vrow(j,i-1),1)*((1-kap)*(vrow(j,i)-vrow(j,i-1))+(1+kap)*(vrow(j,i+1)-vrow(j,i)));

     
      max_abs_u_x = max(abs(ur), abs(ul));
      Fi_x(j,i) = 0.5*(irowr*ur+irowl*ul)- 0.5*max_abs_u_x*(irowr-irowl);
      Fs_x(j,i) = 0.5*(srowr*ur+srowl*ul)- 0.5*max_abs_u_x*(srowr-srowl);
      Fe_x(j,i) = 0.5*(erowr*ur+erowl*ul)- 0.5*max_abs_u_x*(erowr-erowl);
      Fv_x(j,i) = 0.5*(vrowr*ur+vrowl*ul)- 0.5*max_abs_u_x*(vrowr-vrowl);
      
      end

        end
    end

    % Compute numerical fluxes in the y-direction using MUSCL recons.

    for i = 1:Nx+2
        for j = 1:Ny+1

        if any(j<=2) || any(j >=Ny)

        max_abs_v_y = max(abs(v(j,i)), abs(v(j+1,i)));
        Fi_y(j,i) = 0.5*(irow(j,i)*v(j,i)+ irow(j+1,i)*v(j+1,i))- 0.5*max_abs_v_y*(irow(j+1,i) - irow(j,i));  
        Fs_y(j,i) = 0.5*(srow(j,i)*v(j,i)+ srow(j+1,i)*v(j+1,i))- 0.5*max_abs_v_y*(srow(j+1,i) - srow(j,i));  
        Fe_y(j,i) = 0.5*(erow(j,i)*v(j,i)+ erow(j+1,i)*v(j+1,i))- 0.5*max_abs_v_y*(erow(j+1,i) - erow(j,i));  
        Fv_y(j,i) = 0.5*(vrow(j,i)*v(j,i)+ vrow(j+1,i)*v(j+1,i))- 0.5*max_abs_v_y*(vrow(j+1,i) - vrow(j,i));
        else

      vr = v(j+1,i) - 0.25*lim(v(j+2,i)-v(j+1,i),v(j+1,i)-v(j,i),1)*((1-kap)*(v(j+2,i)-v(j+1,i))+(1+kap)*(v(j+1,i)-v(j,i)));
      vl = v(j,i)   + 0.25*lim(v(j+1,i)-v(j,i),v(j,i)-v(j-1,i),1)*((1-kap)*(v(j,i)-v(j-1,i))+(1+kap)*(v(j+1,i)-v(j,i)));

      irowr = irow(j+1,i) -0.25*lim(irow(j+2,i)-irow(j+1,i),irow(j+1,i)-irow(j,i),1)*((1-kap)*(irow(j+2,i)-irow(j+1,i))+(1+kap)*(irow(j+1,i)-irow(j,i)));
      irowl = irow(j,i)   + 0.25*lim(irow(j+1,i)-irow(j,i),irow(j,i)-irow(j-1,i),1)*((1-kap)*(irow(j,i)-irow(j-1,i))+(1+kap)*(irow(j+1,i)-irow(j,i)));

      srowr = srow(j+1,i) -0.25*lim(srow(j+2,i)-srow(j+1,i),srow(j+1,i)-srow(j,i),1)*((1-kap)*(srow(j+2,i)-srow(j+1,i))+(1+kap)*(srow(j+1,i)-srow(j,i)));
      srowl = srow(j,i)  + 0.25*lim(srow(j+1,i)-srow(j,i),srow(j,i)-srow(j-1,i),1)*((1-kap)*(srow(j,i)-srow(j-1,i))+(1+kap)*(srow(j+1,i)-srow(j,i)));

      erowr = erow(j+1,i)  -0.25*lim(erow(j+2,i)-erow(j+1,i),erow(j+1,i)-erow(j,i),1)*((1-kap)*(erow(j+2,i)-erow(j+1,i))+(1+kap)*(erow(j+1,i)-erow(j,i)));
      erowl = erow(j,i)   + 0.25*lim(erow(j+1,i)-erow(j,i),erow(j,i)-erow(j-1,i),1)*((1-kap)*(erow(j,i)-erow(j-1,i))+(1+kap)*(erow(j+1,i)-erow(j,i)));

      vrowr = vrow(j+1,i)  -0.25*lim(vrow(j+2,i)-vrow(j+1,i),vrow(j+1,i)-vrow(j,i),1)*((1-kap)*(vrow(j+2,i)-vrow(j+1,i))+(1+kap)*(vrow(j+1,i)-vrow(j,i)));
      vrowl = vrow(j,i)   + 0.25*lim(vrow(j+1,i)-vrow(j,i),vrow(j,i)-vrow(j-1,i),1)*((1-kap)*(vrow(j,i)-vrow(j-1,i))+(1+kap)*(vrow(j+1,i)-vrow(j,i)));

      max_abs_v_y = max(abs(vr), abs(vl));
      Fi_y(j,i) = 0.5*(irowr*vr+irowl*vl)- 0.5*max_abs_v_y*(irowr-irowl);
      Fs_y(j,i) = 0.5*(srowr*vr+srowl*vl)- 0.5*max_abs_v_y*(srowr-srowl);
      Fe_y(j,i) = 0.5*(erowr*vr+erowl*vl)- 0.5*max_abs_v_y*(erowr-erowl);
      Fv_y(j,i) = 0.5*(vrowr*vr+vrowl*vl)- 0.5*max_abs_v_y*(vrowr-vrowl);

        end
       end
    end

% Update solution

    for i=2:Nx+1
        for j=2:Ny+1

         iirow(j-1,i-1) = irow(j,i) -D_t*(Fi_x(j,i)-Fi_x(j,i-1))/D_x ...
                                    -D_t*(Fi_y(j,i)-Fi_y(j-1,i))/D_y;


         vvrow(j-1,i-1) = vrow(j,i) -D_t*(Fv_x(j,i)-Fv_x(j,i-1))/D_x ...
                                    -D_t*(Fv_y(j,i)-Fv_y(j-1,i))/D_y;

        if  (j == Ny+1) || (i==Nx+1)

        ssrow(j-1,i-1) = srow(j,i) -D_t*(Fs_x(j,i)-Fs_x(j,i-1))/D_x ...
                                   -D_t*(Fs_y(j,i)-Fs_y(j-1,i))/D_y ...
                                   -D_t*b(j-1,i-1)*srow(j,i);

        eerow(j-1,i-1) = erow(j,i) -D_t*(Fe_x(j,i)-Fe_x(j,i-1))/D_x ...
                                   -D_t*(Fe_y(j,i)-Fe_y(j-1,i))/D_y ...
                                   +D_t*b(j-1,i-1)*srow(j,i);
        else

         ssrow(j-1,i-1) = srow(j,i) -D_t*(Fs_x(j,i)-Fs_x(j,i-1))/D_x ...
                                    -D_t*(Fs_y(j,i)-Fs_y(j-1,i))/D_y...
                                    -D_t*b(j,i)*srow(j,i);

        eerow(j-1,i-1) = erow(j,i) -D_t*(Fe_x(j,i)-Fe_x(j,i-1))/D_x ...
                                   -D_t*(Fe_y(j,i)-Fe_y(j-1,i))/D_y ...
                                   +D_t*b(j,i)*srow(j,i);
        end
        end
    end

 %one exit right record exposed peestrian exiting the domain

 for j=NEXY1:NEXY2
    for i=Nx+1:Nx+1
        SUM_EXP_EXIT = SUM_EXP_EXIT + (D_t/D_x)*Fe_x(j,i)*(D_x*D_y);
    end
 end

%RT = iirow + ssrow + eerow + vvrow; 

sumsum=0;
SUMETOTAL=0;


for i=1:Ny
    for j=1:Nx
        sumsum=sumsum+D_x*D_y*(eerow(i,j)+iirow(i,j)+ssrow(i,j)+vvrow(i,j)); %Total pedestrian in the domain

        SUMETOTAL = SUMETOTAL  + D_x*D_y*eerow(i,j); %Exposed pedestrian in the domain

    end
end
fprintf('Mass_Pedestrians = %14.8f\n', sumsum ); 

SUMETOTAL =SUMETOTAL + SUM_EXP_EXIT; %Total exposed pedestrian  over time
fprintf('Mass_Exposed_Pedestrian = %f\n', SUMETOTAL);  
fprintf('***************************\n')

fprintf(fileID,'%8.4f %12.8f %12.8f\n',t*D_t,SUMETOTAL, sumsum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %end of time loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(fileID);
close(writerObj);

figure(2)
surf(x,y,R_new)
title('Total \rho')
colorbar

figure(3)
quiver(x,y,U_new,V_new)
title('Velocity field')

figure(4)
surf(x,y,eerow)
title('Exposed \rho^E')
colorbar

figure(5)
surf(x,y,b)
title('Coef \beta^I')
colorbar

figure(6)
surf(x,y,ssrow)
title('Exposed \rho^S')
colorbar

figure(7)
surf(x,y,iirow)
title('Infected \rho^I')
colorbar
