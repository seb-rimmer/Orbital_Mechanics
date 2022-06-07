%% AE508 - Bang-bang Control Trajectory in LEO
% 
% Design of a VSI trajectory in LEO to get from A to B circular orbits 
% (450 to 500km altitude)
% Earth-based canonical units

%% Errors and Warnings
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

% Canonical Units conversion 
TU = 806.9;             % 1 TU, seconds
DU = 6378;              % Mean radius of earth, km
VU = DU/TU;
mu_earth = 3.98e+14;    % m^3 / s^-2 For ND Units
g_earth  = 9.81;        % Earth g constant, m/s^2

% Choosing spacecraft and propulsion scheme to use
sc_1_mediumsat   = true;
sc_2_cubesat     = false;

ht_prop = false;
lt_prop = true;

% Spacecraft #1 Parameters  
if sc_1_mediumsat == true
% -----------------------------------------------------------------------------
% Roughly based off size and thrust capabilities of a GPS satellite.
%
% Uses the Aerojet-RocketDyne Hydrazine MR-106L for CSI propulsion. 
% https://www.rocket.com/sites/default/files/documents/In-Space%20Data%20Sheets_7.19.21.pdf
%
% Uses NASA Next-C Ion thruster for VSI Propulsion
% https://www1.grc.nasa.gov/wp-content/uploads/NEXT-C_FactSheet_11_1_21_rev4.pdf
% Power limit 0.6 – 7.4 kW
mass  = 500;               % SC mass, kg

if ht_prop == true
isp_s = 230;               % Engine ISP, seconds
T_N   = 20;                % Engine Thrust, N, kg * m/s^2 
end

if lt_prop == true
isp_s = 4220;              % Engine ISP, seconds
T_N   = 0.235;             % Engine Thrust, N, kg * m/s^2 
end

c_ms  = g_earth * isp_s;   % c, m/s
gamma = T_N/mass;
b = gamma*mass/c_ms;
power = 0.5 * T_N * c_ms;
end

% Spacecraft #1 Parameters
if sc_2_cubesat == true
% -----------------------------------------------------------------------------
% Roughly based off size and thrust capabilities of a 12 U CubeSat
%
% Uses the Dawn Aerospace B20 Thruster for CSI propulsion. 
% https://catalog.orbitaltransports.com/b20-thruster-green-propulsion/
%
% Uses L-3 ETI XIPS 13cm Xenon thruster for VSI Propulsion
% https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1459&context=smallsat

mass  = 16;                % SC mass, kg

if ht_prop == true
isp_s = 300;               % Engine ISP, seconds
T_N   = 8;                 % Engine Thrust, N, kg * m/s^2 
end

if lt_prop == true
isp_s = 2350;               % Engine ISP, seconds
T_N   = 0.018;              % Engine Thrust, N, kg * m/s^2 
end

c_ms  = g_earth * isp_s;    % c, m/s
gamma = T_N/mass;
b = gamma*mass/c_ms;
power = 0.5 * T_N * c_ms;
end

% Initial and final altitude conditions 
alt_1 = 400;            % Altitude 1, km
alt_2 = 450;            % Altitude 2, km

% Guest-imate for time of flight using Eq 6.9
t_f_guess = (mu_earth^0.5)/gamma * (((DU+alt_1)*1000)^-0.5 - ((DU+alt_2)*1000)^-0.5 );
t_f_guess_hours = t_f_guess/3600;

% Guest-imate for delta-v using Eq 6.12
dv_tot_guess = gamma * t_f_guess;
dv_tot_guess_kms = dv_tot_guess/1000;


%% Setting up of Initial and Final Conditions

% Non-dimensionalizing SC parameters
T     = T_N * TU^2 / (DU*1000);      
g     = g_earth * TU^2 / (DU*1000);  
isp   = isp_s * 1/TU;                
c     = isp * g;
mu    = 1;

% Initial and final tangential velocities, ND 
v_1   = sqrt(mu_earth/((DU + alt_1)))/1000;
v_2   = sqrt(mu_earth/(DU + alt_2))/1000;

% Initial Conditions 
r_t0     = (DU + alt_1)/DU;                % Radial distance from earth
theta_t0 = 0.0;                            % Angular Position
u_t0     = 0.0;                            % Radial Velocity
v_t0     = sqrt(1/r_t0);                   % Tangential Velocity
% m_t0     = mass;                           % SC Mass
x0       = [r_t0 theta_t0 u_t0 v_t0]; % m_t0];

% Final conditions 
r_tf     = (DU + alt_2)/DU;
u_tf     = 0.0;
v_tf     = sqrt(1/r_tf);
xf       = [r_tf u_tf v_tf];

% Time of flight(s)
t0 = 0;
tf = t_f_guess*2/TU

%% Iterating on initial costate guesses
% Setting ODE and fsolve options
opts_ode = odeset('RelTol',3e-12,'AbsTol',1e-12); % ode tolerances being set
options = optimoptions('fsolve','Display','iter','MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-9,'TolX',1e-10); 

% Determine some random co-state guesses to start iterations with
fuel_burnt = [];

% iterating on tof
while (tf < 250)
    
    % Resetting costate guesses, error, and rho with each iteration
    lam0 = 1e-2 * rand(3, 1)';
    err = 1;
    
    while (err > 1e-2)

        % Setting up params as one object to pass to ode45
        params = [mu T c];

        % Solve OCP
        [lam_new, err] = fsolve(@cost_minU,lam0,options,t0,tf,x0,xf,params,opts_ode);

        lam0 = lam_new;

    end

    % Propagate Min-U Solution
    [t_minU,X_minU] = ode45(@helio_trans,[t0 tf],[x0'; lam0'],opts_ode,params);

    err = 1;            % Resetting error
    tf = tf*1.5;
end

gamma = vecnorm([X_minU(:,6) X_minU(:,7)], 2, 2);       % Units of DU/TU^2
gamma_ms = gamma * TU^2 / (DU*1000);                    % Units of m/s^2
thrust = gamma_ms * mass;
inst_fuel_burnt = [];

for i = 1 : length(gamma_ms)
    b = (mass * gamma_ms(i))^2 / power;
    if i == 1
        dt = t_minU(1);
    else
        dt = t_minU(i) - t_minU(i-1);
    end
    inst_fuel_burnt(end+1) = dt*b;  
end
cum_fuel = cumsum(inst_fuel_burnt);
fuel_burnt = cum_fuel(end);
max_thrust = max(thrust);

% Actual time of flight
actual_tof_TU = t_minU(end);
actual_tof_hours = (actual_tof_TU * TU) / 3600;

%% Plots (final converged solution)
plots(t_minU,X_minU, cum_fuel, thrust);
   
%% Functions
function Xdot = helio_trans(t,X,params)

mu  = params(1);
T   = params(2);
c   = params(3);

r     = X(1);
theta = X(2);
u     = X(3);
v     = X(4);

lam_r     = X(5);
lam_th    = 0;
lam_u     = X(6);
lam_v     = X(7);

% State Differential Equations
xdot = [u;                                  % r_dot
        v/r;                                % th_dot
       (v^2)/r - mu/r^2 - lam_v;            % u_dot
       -u*v/r - lam_u];                     % v_dot
   
lam_dot = [ -(-lam_th*v/r^2 - lam_u*(v^2)/(r^2) + 2*lam_u*mu/r^3 + lam_v*u*v/r^2);
            -(lam_r - lam_v*v/r);
            -(lam_th/r + 2*lam_u*v/r - lam_v*u/r);
            ];
       
Xdot = [xdot; lam_dot;];

end

function err = cost_minU(lam0_guess,t0,tf,x0,xf,params,opts_ode)

[t,X] = ode45(@helio_trans,[t0 tf],[x0'; lam0_guess'], opts_ode, params);

err = [X(end,1) 0 X(end,3) X(end,4)]' - [xf(1) 0 0 sqrt(1/xf(1))]';

end

function plots(t_minU,X_minU, cum_fuel, thrust)

%% Trajectory
figure(2)
radial_distance = X_minU(:,1);
angular_position = X_minU(:,2);
polarplot(angular_position, radial_distance)
pax = gca;
pax.RLim = [X_minU(1,1)*0.98, X_minU(end,1)*1.01];
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

% subplot 235
% plot(t_minU,X_minU(:,5),'b-','Linewidth',2)
% grid minor
% ylabel('Mass (kg)')
% xlabel('Time (TU)')
% title('State x_5 - Spacecraft mass')

%% Co-States
figure(4)
subplot 311
plot(t_minU,X_minU(:,5),'b-','Linewidth',2)
grid minor
ylabel('\lambda_r')
title('Radial Distance Costate \lambda_r')

subplot 312
plot(t_minU,X_minU(:,6),'b-','Linewidth',2)
grid minor
ylabel('\lambda_u')
title('Radial Velocity Costate \lambda_u')

subplot 313
plot(t_minU,X_minU(:,7),'b-','Linewidth',2)
grid minor
ylabel('\lambda_v')
xlabel('Time (TU)')
title('Tangential Velocity Costate \lambda_v')

%% Control - Can plot from the co-states
figure(5)
title('Control Vectors for final solution')

subplot 411
plot(t_minU, -X_minU(:,7),'b-','Linewidth',2)
ylabel('$u_r (AU/TU^2)$')
xlabel('Time (TU)')
grid minor
title('Control - Radial Thrust Component')

subplot 412
plot(t_minU, -X_minU(:,6),'b-','Linewidth',2)
ylabel('$u_\theta (AU/TU^2)$')
xlabel('Time (TU)')
grid minor
title('Control - Angular Thrust Component')

subplot 413
plot(t_minU, thrust,'b-','Linewidth',2)
ylabel('Thrust - N')
xlabel('Time (TU)')
grid minor
title('Required Thrust Profile')

subplot 414
plot(t_minU, cum_fuel,'b-','Linewidth',2)
ylabel('mass (kg)')
xlabel('Time (TU)')
grid minor
title('Cummulative Expended Fuel')
end
