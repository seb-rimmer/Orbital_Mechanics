%% AE508 Final Project
% -------------------------------------------------------------------------
% Comparison of VSI vs CSI propulsion schemes in LEO to return to a nominal
% trajectory.
%
% Modelling a simple orbit raise from A to B. Fixed final time problem.
%
% Earth-based canonical units used
% -------------------------------------------------------------------------
%% System Variables
close all; clc;

% Canonical Units conversion 
M_e = 6.674e-11;                    % mass of Earth, kg
G   = 5.972e+24;                    % Gravitational parameter, N m^2 kg^-2
mu_earth = G * M_e;                 % m^3 / s^-2 For ND Units

DU  = 6378;                         % Mean radius of earth, km
TU = sqrt((DU*1000)^3/(mu_earth));  % Time unit from canonical units equation
VU = DU/TU;

g_earth  = 9.81;                    % Earth g constant, m/s^2

%% Choosing spacecraft and propulsion scheme to use. 
% We solve for one of 4 possible scenarios:
% 
%   1. Medium spacecraft, CSI high-thrust chemical propulsion
%   2. CubeSat, CSI high-thrust chemical propulsion
%
%   3. Medium spacecraft, VSI low-thrust electric propulsion
%   4. CubeSat, VSI low-thrust electric propulsion
%
scenario = 1;

if mod(scenario, 2)
    sc_1_mediumsat = true;
    sc_2_cubesat   = false;
else
    sc_1_mediumsat = false;
    sc_2_cubesat   = true; 
end
if scenario <= 2
    ht_prop = true;
    lt_prop = false;
    lam_guess_mag = 1e-1;
else
    ht_prop = false;
    lt_prop = true;
    lam_guess_mag = 1e-2;
end

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

c_ms  = g_earth * isp_s;   % c, m/s
power = 0.5 * T_N * c_ms;  % Resulting propulsion system power requirement
b = gamma*mass/c_ms;
end

gamma = T_N/mass;

end
% Spacecraft #2 Parameters
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
power = 0.5 * T_N * c_ms;  % Resulting propulsion system power requirement
end

c_ms  = g_earth * isp_s;    % c, m/s
gamma = T_N/mass;
b = gamma*mass/c_ms;
end

%% Initial and final conditions 
alt_1 = 400;            % Altitude 1, km
alt_2 = 450;            % Altitude 2, km

% Guest-imate for time of flight using Eq 6.9
t_f_guess = (mu_earth^0.5)/gamma * (((DU+alt_1)*1000)^-0.5 - ((DU+alt_2)*1000)^-0.5 );
t_f_guess_hours = t_f_guess/3600;

% Guest-imate for delta-v using Eq 6.12
dv_tot_guess = gamma * t_f_guess;
dv_tot_guess_kms = dv_tot_guess/1000;

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
m_t0     = mass;                           % SC Mass
x0       = [r_t0 theta_t0 u_t0 v_t0];
if ht_prop
   x0(end+1) = m_t0; 
end

% Final conditions 
r_tf     = (DU + alt_2)/DU;
u_tf     = 0.0;
v_tf     = sqrt(1/r_tf);
xf       = [r_tf u_tf v_tf];

% Propagator limits
tf_guess_enlarge_factor = [3 10 2.3 2];
tf_upper_limit          = [4.5 4 210 70]; 
tf_growth_factor        = [1.05 1.2 1.02 1.02];

tf_guess_enlarge_factor = tf_guess_enlarge_factor(scenario);
tf_upper_limit          = tf_upper_limit(scenario);
tf_growth_factor        = tf_growth_factor(scenario);

% Time of flight(s)
t0 = 0;
tf = t_f_guess*tf_guess_enlarge_factor/TU;

%% Iterating on initial costate guesses
% Setting ODE and fsolve options
opts_ode = odeset('RelTol',3e-12,'AbsTol',1e-12); % ode tolerances being set
options = optimoptions('fsolve','Display','iter','MaxFunEvals',100000,'MaxIter',100000,'TolFun',1e-8,'TolX',1e-10); 

% Determine some random co-state guesses to start iterations with
lam0 = lam_guess_mag * rand(3, 1)';
err = 1;
rho = 1;
fuel_burnt = [];

% iterating on tof
while (tf < tf_upper_limit)
    
    % Resetting error, and rho with each iteration
    lam0 = lam_guess_mag * rand(3, 1)';
    err = 1;
    rho = 1;
    
    while (rho > 1e-4)

        while (err > 1e-2)

            % Setting up params as one object to pass to ode45
            params = [mu T rho c ht_prop];

            % Solve OCP
            [lam_new, err] = fsolve(@cost_minU,lam0,options,t0,tf,x0,xf,params,opts_ode);

            lam0 = lam_new;

        end

        % Propagate Min-U Solution
        [t_minU,X_minU] = ode45(@leo_transfer,[t0 tf],[x0'; lam0'],opts_ode,params);

        rho = rho * 0.3;    % Sweeping rho down for next iteration
        err = 1;            % Resetting error
    end
       
    if ht_prop == true
        
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
        rho_str = sprintf('Time of flight = %.3f TU',tf);

        figure(1);
        plot(t_minU, s_t, 'DisplayName', rho_str);    
        xlabel('Time (TU)','FontSize', 15)
        xlim([0 t_minU(end)])
        ylabel('Switch Function Magnitude (Non-Di)','FontSize', 15)
        title(sprintf(fig_title, scenario),'FontSize', 18)       
%         legend()
        grid on
        grid minor
        hold all

        % Plotting new thrust profile with each iteration
        figure(2);
        plot(t_minU, thrust, 'DisplayName', rho_str);
        xlabel('Time (TU)','FontSize', 15)
        xlim([0 t_minU(end)])
        ylim([0 25])
        ylabel('Thrust (N)','FontSize', 15)
        fig_title = 'Thurst Profile for CSI Scenario %d'; 
        title(sprintf(fig_title, scenario),'FontSize', 18)
%         legend()
        grid on
        grid minor
        hold all

        fuel_mass = mass - X_minU(end,5);
        fuel_burnt(end+1, 1) = t_minU(end);
        fuel_burnt(end, 2) = fuel_mass;
        cum_fuel = nan;
    else
        gamma = vecnorm([X_minU(:,6) X_minU(:,7)], 2, 2);   % Units of DU/TU^2
        gamma_ms = gamma * TU^2 / (DU*1000);                % Units of m/s^2
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
        max_thrust = max(thrust);
        
        fuel_burnt(end+1, 1) = t_minU(end);
        fuel_burnt(end, 2) = cum_fuel(end);
        
        u_1 = nan;
        u_2 = nan;
    end
    
    tf = tf*tf_growth_factor;
end

if ht_prop == true
    file_name = 'scenario_%d_results/scenario_%d_switch_function.png';
    output_file_name_sf = sprintf(file_name, scenario, scenario);
    saveas(figure(1), output_file_name_sf)
    
    file_name = 'scenario_%d_results/scenario_%d_thrust_profile.png';  
    output_file_name_tp = sprintf(file_name, scenario, scenario);
    saveas(figure(2), output_file_name_tp)
end
  
figure(6);
scatter(fuel_burnt(:,1)',fuel_burnt(:,2),20,'k','x');
grid minor
ylabel('Cost (Fuel mass, kg)', 'FontSize', 15)
xlabel('Time of flight (TU)', 'FontSize', 15)
fig_title = 'Cost vs. Time of Flight for Scenario %d'; 
title(sprintf(fig_title, scenario), 'FontSize', 18)
file_name = 'scenario_%d_results/scenario_%d_cost_vs_tof.png';
output_file_name = sprintf(file_name, scenario, scenario);
saveas(figure(6),output_file_name);
    


%% Plots
plot_trajectory(X_minU, scenario);
plot_states(t_minU,X_minU, ht_prop);
% plot_costates(t_minU,X_minU);
plot_controls(t_minU,X_minU, u_1, u_2, ht_prop, thrust, cum_fuel)

%% Functions
function d = delta(rho, lam_u, lam_v)

p = -[lam_u, lam_v];
S = vecnorm(p, 2, 2) - 1; 

d = 0.5 * (1 + tanh(S/rho));

end

function Xdot = leo_transfer(t,X,params)

mu  = params(1);
T   = params(2);
rho = params(3);
c   = params(4);
ht_prop = params(5);

r     = X(1);
theta = X(2);
u     = X(3);
v     = X(4);

if ht_prop 
    m     = X(5);
    lam_r     = X(6);
    lam_th    = 0;
    lam_u     = X(7);
    lam_v     = X(8);
    
    u_1 = - lam_u./vecnorm([lam_u lam_v], 2, 2);
    u_2 = - lam_v./vecnorm([lam_u lam_v], 2, 2);
    d   =   delta(rho, lam_u, lam_v);

else
    lam_r     = X(5);
    lam_th    = 0;
    lam_u     = X(6);
    lam_v     = X(7);
end

% State Differential Equations
if ht_prop == true
    xdot = [u;
            v/r;
           (v^2)/r - mu/r^2 + T/m * d * u_1;
           -u*v/r + T/m * d * u_2;
           -(T/c) * d;];     % Added state equation to describe mass
else
    xdot = [u;                                  % r_dot
            v/r;                                % th_dot
           (v^2)/r - mu/r^2 - lam_v;            % u_dot
           -u*v/r - lam_u];                     % v_dot
end

lam_dot = [ -(-lam_th*v/r^2 - lam_u*(v^2)/(r^2) + 2*lam_u*mu/r^3 + lam_v*u*v/r^2);
            -(lam_r - lam_v*v/r);
            -(lam_th/r + 2*lam_u*v/r - lam_v*u/r);
            ];
       
Xdot = [xdot; lam_dot;];

end

function err = cost_minU(lam0_guess,t0,tf,x0,xf,params,opts_ode)

[t,X] = ode45(@leo_transfer,[t0 tf],[x0'; lam0_guess'], opts_ode, params);

err = [X(end,1) 0 X(end,3) X(end,4)]' - [xf(1) 0 0 sqrt(1/xf(1))]';

end

function plot_trajectory(X_minU, scenario)
    figure(7)
    radial_distance = X_minU(:,1);
    angular_position = X_minU(:,2);
    polarplot(angular_position, radial_distance)
    pax = gca;
    pax.RLim = [X_minU(1,1)*0.98, X_minU(end,1)*1.01];
    
    fig_title = 'Polar Plot of Trajectory for Scenario %d'; 
    title(sprintf(fig_title, scenario), 'FontSize', 18)
    
    file_name = 'scenario_%d_results/scenario_%d_trajectory.png';
    output_file_name = sprintf(file_name, scenario, scenario);
    saveas(figure(7),output_file_name);

end

function plot_states(t_minU,X_minU,ht_prop)
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

    if ht_prop
        subplot 235
        plot(t_minU,X_minU(:,5),'b-','Linewidth',2)
        grid minor
        ylabel('Mass (kg)')
        xlabel('Time (TU)')
        title('State x_5 - Spacecraft mass')
    end
end

function plot_costates(t_minU,X_minU)
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
end

function plot_controls(t_minU,X_minU, u_1, u_2, ht_prop, thrust, cum_fuel)

%% Control - Can plot from the co-states
figure(5)
title('Control Vectors for final solution')

if ht_prop
    subplot 211
    plot(t_minU, u_1,'b-','Linewidth',2)
    ylabel('$u_r (AU/TU^2)$')
    xlabel('Time (TU)')
    grid minor
    title('Control - Radial Thrust Component')

    subplot 212
    plot(t_minU, u_2,'b-','Linewidth',2)
    ylabel('$u_\theta (AU/TU^2)$')
    xlabel('Time (TU)')
    grid minor
    title('Control - Angular Thrust Component')
else

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

end
