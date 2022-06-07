%% AE508 - Bang-bang Control Trajectory to Mars
% 
% You are required to design a bang-bang control to send a spacecraft from 
% Earth to Mars with a time of flight of 300 days (1 TU = 58.13 days). Assume
% Earth and Mars orbit the Sun in co-planar circular orbits of radii 1 AU and 
% 1.524 AU respectively. 

%% Errors
%{
Warning: Trust-region-dogleg algorithm of FSOLVE cannot handle non-square systems; using Levenberg-Marquardt algorithm instead. 
> In fsolve (line 316)
  In AE508_HW4_Rimmer (line 60) 

                                        First-Order                    Norm of 
 Iteration  Func-count    Residual       optimality      Lambda           step
     0           5        0.678204            14.3         0.01
     1          20        0.678204            6.17        1e+08    1.62021e-07

No solution found.

fsolve stopped because the relative size of the current step is less than the
value of the step size tolerance, but the vector of function values
is not near zero as measured by the value of the function tolerance. 
  
- Is this a canonical units problem? Unsure on conversion on line 60-63
- Have we set up rho function correctly? p correct?
- Lines 105/107, not including m or theta yet, not sure if ok or not
- Is error done right?
%}

%% System Variables
close all; clc;

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
x0 = [r_t0 theta_t0 u_t0 v_t0 m_t0];

% Final conditions 
r_tf     = 1.524;
u_tf     = 0.0;
v_tf     = sqrt(1/r_tf);
xf       = [r_tf nan v_tf u_tf nan];

% Time of flight(s)
t0 = 0;
tf = 300 / TU;

% Need to non-dimensionalize thrust parameters
T_N   = 0.235;          % Thrust, N, kg * m/s^2 
g_E   = 9.81;           % Earth G constant, m/s^2
isp_s = 4155;           % ISP, s
c_ms  = g_E * isp_s;    % c, m/s

T     = T_N * (TU * 24*3600)^2 / (DU);      % DU should be in m here?
g     = g_E * (TU * 24*3600)^2 / (DU);      % Does not give T = 39.3123 or
isp   = isp_s * 1/(TU * 24*3600);             % c = 1.3644
c     = isp * g;

% Determine Co-state guesses to work with at different orders of magnitude
lam0_guess = [-2; 0.5; -2;];

%% Solving

% Setting ODE and fsolve options
opts_ode = odeset('RelTol',3e-14,'AbsTol',1e-14); % ode tolerances being set
options = optimoptions('fsolve','Display','iter','MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-8,'TolX',1e-10); 

rho_shrink_factor = 0.3;

% Setting up figure for Switch Function and Thrust Profile plots
figure(1)

while (rho > 1e-3)
    
    % Setting up params (as one object to pass to ode45)
    params = [mu T rho c];
    
    % Solve OCP
    [lam0,fval] = fsolve(@cost_minU,lam0_guess,options,t0,tf,x0,xf,params,opts_ode);

    % Propagate Min-U Solution
    [t_minU,X_minU] = ode45(@helio_trans,[t0 tf],[x0'; lam0],opts_ode,params);
     
    % Working out controls from co-states
    lam_u = X_minU(:,7);
    lam_v = X_minU(:,8);

    % Needs to be vec norm
    u_1 = -(lam_u./vecnorm([lam_u lam_v], 2, 2));
    u_2 = -(lam_v./vecnorm([lam_u lam_v], 2, 2));
    
    % Determining switch function and thrust profile for iteration
    p = -[lam_u, lam_v];
    s_t = vecnorm(p, 2, 2) - 1; 
    throttle = arrayfun(@delta,  rho* ones(length(lam_u),1), lam_u, lam_v);
    thrust = T_N .* throttle;
    
    % Plotting new switch function profile with each iteration
    rho_str = sprintf('\\rho = %.3f',rho);
    
    subplot 211
    plot(t_minU, s_t, 'DisplayName', rho_str);    
    xlabel('Time')
    xlim([0 t_minU(end)])
    ylabel('Switch Function\n Magnitude (Non-Di)')
    title('Switch Function')
    legend()
    grid on
    grid minor
    hold all

    % Plotting new thrust profile with each iteration
    subplot 212
    plot(t_minU, thrust, 'DisplayName', rho_str);
    xlabel('Time')
    xlim([0 t_minU(end)])
    ylabel('Thrust (N)')
    title('Thurst Profile')
    legend()
    grid on
    grid minor
    hold all

    % Sweeping rho down for next iteration
    rho = rho * rho_shrink_factor;
end

%% Other Plots (final converged solution)
plots(t_minU,X_minU, u_1, u_2)

%% Functions
function d = delta(rho, lam_u, lam_v)

S = norm([lam_u lam_v]) - 1;

d = 0.5 * (1 + tanh(S/rho));

end

function Xdot = helio_trans(t,X,params)

mu  = params(1);
T   = params(2);
rho = params(3);
c   = params(4);

r     = X(1);
theta = X(2);
u     = X(3);
v     = X(4);
m     = X(5);

lam_r     = X(6);
lam_th    = 0; %X(7);
lam_u     = X(7);
lam_v     = X(8);
% lam_m     = X(10);    % Can comment this out along with line 190

u_1 = - lam_u./vecnorm([lam_u lam_v], 2, 2);
u_2 = - lam_v./vecnorm([lam_u lam_v], 2, 2);
d   =   delta(rho, lam_u, lam_v);

% State Differential Equations
xdot = [u;
        v/r;
       (v^2)/r - mu/r^2 + T/m * d * u_1;
       -u*v/r + T/m * d * u_2;
       -(T/c) * d;];     % Added state equation to describe mass
   
% lam_th = 0;
lam_dot = [ -(-lam_th*v/r^2 - lam_u*(v^2)/(r^2) + 2*lam_u*mu/r^3 + lam_v*u*v/r^2);
            
            % We know lam_theta is zero, so can be removed
%             0;
            -(lam_r - lam_v*v/r);
            -(lam_th/r + 2*lam_u*v/r - lam_v*u/r);
            
            % Should be able to comment out lam_m costate, doesn't appear 
            % in other costate equations.
%             -(-T/m^2 * d - lam_u*T/m^2 * d * u_1 - lam_v*T/m^2 * d * u_2);
            ];
       
Xdot = [xdot; lam_dot;];

end

function err = cost_minU(lam0_guess,t0,tf,x0,xf,params,opts_ode)

[t,X] = ode45(@helio_trans,[t0 tf],[x0'; lam0_guess], opts_ode, params);

err = [X(end,1) 0 X(end,3) X(end,4)]' - [xf(1) 0 0 sqrt(1/xf(1))]';

end

function plots(t_minU,X_minU, u_1, u_2)

%% Trajectory
figure(2)
radial_distance = X_minU(:,1);
angular_position = X_minU(:,2);
polarplot(angular_position, radial_distance)
pax = gca;
pax.RLim = [0.8, 1.6];
title('Polar Plot of Trajectory')

%% States
figure(3)
subplot 231
plot(t_minU,X_minU(:,1),'b-','Linewidth',2)
grid minor
ylabel('Radius (AU)')
xlabel('Time (TU)')
title('State x_1 - Radial Distance')

subplot 232
plot(t_minU,X_minU(:,2),'b-','Linewidth',2)
grid minor
ylabel('\\theta (Radians)')
xlabel('Time (TU)')
title('State x_2 - Angular Position')

subplot 233
plot(t_minU,X_minU(:,3),'b-','Linewidth',2)
grid minor
ylabel('Radial Velocity (Radians/TU)')
xlabel('Time (TU)')
title('State x_3 - Radial Velocity')

subplot 234
plot(t_minU,X_minU(:,4),'b-','Linewidth',2)
grid minor
ylabel('Velocity (AU/TU)')
xlabel('Time (TU)')
title('State x_4 - Tangential Velocity')

subplot 235
plot(t_minU,X_minU(:,5),'b-','Linewidth',2)
grid minor
ylabel('Mass (kg)')
xlabel('Time (TU)')
title('State x_5 - Spacecraft mass')

%% Co-States
figure(4)
subplot 311
plot(t_minU,X_minU(:,6),'b-','Linewidth',2)
grid minor
ylabel('\lambda_r')
title('Radial Distance Costate \lambda_r')

subplot 312
plot(t_minU,X_minU(:,7),'b-','Linewidth',2)
grid minor
ylabel('\lambda_u')
title('Radial Velocity Costate \lambda_u')

subplot 313
plot(t_minU,X_minU(:,8),'b-','Linewidth',2)
grid minor
ylabel('\lambda_v')
xlabel('Time (TU)')
title('Tangential Velocity Costate \lambda_v')

%% Control - Can plot from the co-states
figure(5)
title('Control Vectors for final solution')

subplot 211
plot(t_minU, u_1,'b-','Linewidth',2)
ylabel('$u_r (AU/TU^2)$','Interpreter','latex')
xlabel('Time (TU)')
grid minor
title('Control - Radial Thrust Component')

subplot 212
plot(t_minU, u_2,'b-','Linewidth',2)
ylabel('$u_\theta (AU/TU^2)$','Interpreter','latex')
xlabel('Time (TU)')
grid minor
title('Control - Angular Thrust Component')

end
