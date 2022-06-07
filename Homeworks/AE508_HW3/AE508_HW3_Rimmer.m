%% AE508 HW 3 - Heliocentric Transfer from Earth to Mars
% 
%  You are required to design a spacecraft transfer trajectory from Earth to
%  Mars with a time-of-flight of 150 days (1 TU = 58.13 days). Assume that 
%  Earth and Mars orbit the Sun in co-planar, circular orbits, of radius 
%  1 AU and 1.5 24 AU respectively.
% 
%  This is a minimum time problem, corresponding to   

clear all; close all; clc;

%% System Variables

% System Constansts
TU = 58.13;   % 1 TU, days
mu = 1;       % For ND Units

% Initial conditions 
r_t0     = 1.0;         % Radial distance from sun
v_t0     = 1.0;         % Tangential Velocity
theta_t0 = 0.0;         % Angular Position
u_t0     = 0.0;         % Radial Velocity
x0 = [r_t0 v_t0 theta_t0 u_t0];

% Final conditions 
r_tf     = 1.524;
v_tf     = sqrt(1/r_tf);
u_tf     = 0.0;
xf = [r_tf v_tf u_tf];

% Time of flight(s)
time_of_flight = [70, 100, 150, 200, 360, 600];
t0 = 0;
tf = time_of_flight(3) / TU;

% Sets of Co-state guesses to work with at different orders of magnitude
% lam0_guess = [-0.1; -0.1; -0.1; -0.1; -0.1; -0.1];
lam0_guess = [-0.01; -0.01; -0.01; -0.01;]; %; -0.01; -0.01];
% lam0_guess = [0.01; 0.01; 0.01; 0.01; 0.01; 0.01];

%% Solving

% Setting ODE and fsolve options
opts_ode = odeset('RelTol',3e-14,'AbsTol',1e-14); % ode tolerances being set
options = optimoptions('fsolve','Display','iter','MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8,'TolX',1e-10); % fsolve options

[lam0,fval] = fsolve(@cost_minU,lam0_guess,options,t0,tf,x0,xf,mu,opts_ode);

% Propagate Min-U Solution
[t_minU,X_minU] = ode45(@helio_trans,[t0 tf],[x0'; lam0],opts_ode,mu);

% Compute Hamiltonian
H = hamiltonian(t_minU,X_minU,mu);

% Computing Cost
L = accum_cost(X_minU);

% Accumulative Cost

%% Plots

% States, Costates(lambdas), Controls, Accumultive Cost
plots(t_minU,X_minU,H,L)

%% Functions
function Xdot = helio_trans(t,X,mu)

x   = X(1:4);   % of the form  x_state = [r_t0 v_t0 theta_t0 u_t0];
lam = X(5:8);

% State Differential Equations
xdot = [x(4);
        -x(4)*x(2)/x(1) - lam(2);
        x(2)/x(1);
        (x(2)^2)/(x(1)) - mu/x(1)^2 - lam(4)];
    
% Costate Differential Equations
lam_dot = [ -(lam(2)*x(4)*x(2)/x(1)^2 + lam(3)*x(2)/x(1)^2 - lam(4)*(x(2)^2)/(x(1)^2) + lam(4)*2*mu/(x(1)^3));
            -(-lam(2)*x(4)/x(1) + lam(3)/x(1) + 2*lam(4)*x(2)/x(1));
            0;
            -(lam(1) - lam(1)*x(2)/x(1));];
         
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

function err = cost_minU(lam0_guess,t0, tf,x0,xf,mu,opts_ode)

[t,X] = ode45(@helio_trans,[t0 tf],[x0'; lam0_guess], opts_ode, mu);

H = hamiltonian(t,X,mu);

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
