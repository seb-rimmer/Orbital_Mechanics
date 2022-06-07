%% AE508 HW 3 - Heliocentric Transfer from Earth to Mars
% 
%  You are required to design a spacecraft transfer trajectory from Earth to
%  Mars with a time-of-flight of 150 days (1 TU = 58.13 days). Assume that 
%  Earth and Mars orbit the Sun in co-planar, circular orbits, of radius 
%  1 AU and 1.524 AU respectively.
% 
%  This is a minimum time problem, corresponding to   

clear all; close all; clc;

%% System Variables

% System Constansts
TU = 58.13;             % 1 TU, days
DU = 149.6e+9;          % 1 AU, km, Sun-Earth Distance
mu = 1;                 % For ND Units
rho = 1;

% Initial conditions 
r_t0     = 1.0;         % Radial distance from sun
theta_t0 = 0.0;         % Angular Position
u_t0     = 0.0;         % Radial Velocity
v_t0     = 1.0;         % Tangential Velocity
m_t0     = 500;         % SC Mass
x0 = [r_t0 v_t0 theta_t0 u_t0];

% Final conditions 
r_tf     = 1.524;
u_tf     = 0.0;
v_tf     = sqrt(1/r_tf);
xf = [r_tf v_tf u_tf];

% Time of flight(s)
tf = 300 / TU;

% Need to non-dimensionalize thrust parameters
T_N   = 0.235;          % Thrust, N, kg * m/s^2 
g_E   = 9.81;           % Earth G constant, m/s^2
isp_s = 4155;           % ISP, m/s

T     = T_N * (TU * 24*3600)^2 / (DU);      % DU should be in m here?
g     = g_E * (TU * 24*3600)^2 / (DU);      % Does not give T = 39.3123 or
isp   = isp_s / (TU * 24*3600);             % c = 1.3644
c     = isp * g;

% Determine Co-state guesses to work with at different orders of magnitude
lam0_guess = [-2; -0.1; -0.5; -2;];

%% Solving

% Setting ODE and fsolve options
opts_ode = odeset('RelTol',3e-14,'AbsTol',1e-14); % ode tolerances being set
options = optimoptions('fsolve','Display','iter','MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8,'TolX',1e-10); % fsolve options

[lam0,fval] = fsolve(@cost_minU,lam0_guess,options,0,tf,x0,xf,mu,T,c,opts_ode);

% Propagate Min-U Solution
[t_minU,X_minU] = ode45(@helio_trans,[0 tf],[x0'; lam0],opts_ode,mu,T,c);

% Compute Hamiltonian
% H = hamiltonian(t_minU,X_minU,mu);

% Computing Cost
% L = accum_cost(X_minU);

% Accumulative Cost

%% Plots

% States, Costates(lambdas), Controls, Accumultive Cost
% plots(t_minU,X_minU,H,L)

%% Functions
function d = delta_t(rho, lam_u, lam_v)

S = norm(-[lam_u lam_v]) - 1;

d = 0.5 * (1 + tanh(S/rho));

end

function Xdot = helio_trans(t,X,mu,T,c)

r     = X(1);
theta = X(2);
u     = X(3);
v     = X(4);
m     = 500;

lam_r     = X(5);
lam_th    = X(6);
lam_u     = X(7);
lam_v     = X(8);
% lam_m     = X(9);

u_1 = - lam_u/norm([lam_u lam_v]);
u_2 = - lam_v/norm([lam_u lam_v]);
d   =   delta_t(1, lam_u, lam_v);

% State Differential Equations
xdot = [u;
        v/r;
       (v^2)/r - mu/r^2 + T/m * d * u_1;
       -u*v/r + T/m * d * u_2;];
%        -(T/c) * d;];     % Added state equation to describe mass

lam_dot = [ -(-lam_th*v/r^2 - lam_u*(v^2)/(r^2) + 2*lam_u*mu/r^3 + lam_v*u*v/r^2);
            0;
            -(lam_r - lam_v*v/r);
            -(lam_th/r + 2*lam_u*v/r - lam_v*u/r);];
%             -(-T/m^2 * d - lam_u*T/m^2 * d * u_1 - lam_v*T/m^2 * d * u_2);];

         
Xdot = [xdot; lam_dot;];

end

function L = accum_cost(X_minU)

L = zeros(length(X_minU), 1);

for i = 1:length(X_minU)
    
    % 
    val_to_add = 0.5*(X_minU(6)^2 + X_minU(8)^2);
    
    % Hamiltonian - is wrong, fix lambda
    L(i) = 0.5*(X_minU(i, 6)^2 + X_minU(i, 8)^2);
    
end

L = cumsum(L);

end

function err = cost_minU(lam0_guess,t0,tf,x0,xf,mu,T,c,opts_ode)

[t,X] = ode45(@helio_trans,[0 tf],[x0'; lam0_guess], opts_ode, mu,T,c);

% H = hamiltonian(t,X,mu);

% Our state function error is only between r, v, and u, NOT theta, hence
% 3x1 sized vector
err = [X(end,1:2) X(end, 4)]' - xf';

end

function H = hamiltonian(t,X,mu)

x   = X(:,1:4);
lam = X(:,5:8);

for i = 1:length(X)
    
    % State Differential Equations
    xdot = [x(4);
            -x(4)*x(2)/x(1) - lam(2);
            x(2)/x(1);
            (x(2)^2)/(x(1)) - mu/x(1)^2 - lam(4)];

    % Hamiltonian - is wrong, fix lambda
    H(i) = 0.5*(lam(2)^2 + lam(4)^2) + lam(i,:)*xdot;
    
end

end

function plots(t_minU,X_minU,H,L)

%% Trajectory
figure(1)
radial_distance = X_minU(:,1);
angular_position = X_minU(:,3);
polarplot(angular_position, radial_distance)
pax = gca;
pax.RLim = [0.8, 1.6];

%% States
figure(2)
subplot 221
plot(t_minU,X_minU(:,1),'b-','Linewidth',2)
ylabel('x_1')
xlabel('Time')
title('State x_1 - Radial Distance')

subplot 222
plot(t_minU,X_minU(:,2),'b-','Linewidth',2)
ylabel('x_2')
xlabel('Time')
title('State x_2 - Tangential Velocity')

subplot 223
plot(t_minU,X_minU(:,3),'b-','Linewidth',2)
ylabel('x_3')
xlabel('Time')
title('State x_3 - Angular Position')

subplot 224
plot(t_minU,X_minU(:,4),'b-','Linewidth',2)
ylabel('x_4')
xlabel('Time')
title('State x_4 - Radial Velocity')

%% Co-States
figure(3)

subplot 231
plot(t_minU,X_minU(:,5),'b-','Linewidth',2)
ylabel('\lambda_1')
xlabel('time')
subplot 232
plot(t_minU,X_minU(:,6),'b-','Linewidth',2)
ylabel('\lambda_2')
xlabel('time')
title('Costates')
subplot 233
plot(t_minU,X_minU(:,7),'b-','Linewidth',2)
ylabel('\lambda_3')
xlabel('time')
subplot 234
plot(t_minU,X_minU(:,8),'b-','Linewidth',2)
ylabel('\lambda_4')
xlabel('time')
% subplot 235
% plot(t_minU,X_minU(:,9),'b-','Linewidth',2)
% ylabel('\lambda_5')
% xlabel('time')
% subplot 236
% plot(t_minU,X_minU(:,10),'b-','Linewidth',2)
% ylabel('\lambda_6')
% xlabel('time')

%% Control - Can plot from the co-states
figure(4)

subplot 211
plot(t_minU, -X_minU(:,8),'b-','Linewidth',2)
ylabel('u_r (AU/TU^2)')
xlabel('Time')
title('Control - Radial Thrust Component')

subplot 212
plot(t_minU, -X_minU(:,6),'b-','Linewidth',2)
ylabel('u_\theta (AU/TU^2)')
xlabel('Time')
title('Control - Angular Thrust Component')

% Accumulative Cost
figure(5)
subplot 111
plot(t_minU, L,'b-','Linewidth',2)
ylabel('Accumulative cost (AU/TU^2)')
xlabel('Time')
title('Accumulative Cost over time')

end

%% HW3 Helper Code
