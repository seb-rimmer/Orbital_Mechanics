%% AE 508 - Lecture 7
% Linear Oscillator - Min Time & Min Control (Fuel)

%% Notes
%  -----------------------------------------------------------------------
%   We are propagating 4 differential equations, seeing if it hits the
%   final target (ie our desired state). Not really a physcial
%   interpretation of the co-states; adjoint variable that helps us solve
%   the problem.  
% 
%   Remember, the lagrange multiplier doesn't have a physical meaning but
%   it can provide a "measure of sensitivity to small changes in initial 
%   conditions"
%   Kind of  need to make a 'random guess' at first to propagate problem
%   with - and the more complicated the problem, the harder it is to define
%   what a "good" random guess looks like. 



%   Can compare performance of different integrator regimes and compare
%   difference in the Hamiltonian error output. 




%%
clear all; close all; clc;


minT = true;
minU = true;
analytical = false;

opts_ode = odeset('RelTol',1e-12,'AbsTol',1e-14); % ode tolerances being set
options = optimoptions('fsolve','Display','iter','MaxFunEvals',100,'MaxIter',100,'TolFun',1e-10,'TolX',1e-14); % fsolve options

% Initialize
t_minU = [];
t_minT = [];
X_minU = [];
X_minT = [];
H_minU = [];
H_minT = [];

% Constants, like mu in orbital mechanics problem to describe dynamics
a = 0.5;
b = 0.1;

% Initial Conditions
x0 = [1 1];
xf = [2 3];

% Time
t0 = 0;
tf = 12;

%% Numerical Solution (Min u)

if minU == true
    
    lam0_guess = [0.1; 0.1];
    
    [lam0,fval] = fsolve(@cost_minU,lam0_guess,options,t0,tf,x0,xf,a,b,opts_ode);
    
    % Propagate Min-U Solution
    [t_minU,X_minU] = ode45(@linosc,[t0 tf],[x0'; lam0],opts_ode,a,b);
    
    % Compute Hamiltonian
    H = hamiltonian(t_minU,X_minU,a,b);
    H_minU = H' + 0.5.*(-X_minU(:,4)).^2;
    
end

%% Numerical Solution (Min Time)

if minT == true
    
    lam0_tf_guess = [0.1; 0.1; 10];   % We are making a guess of 10 for time (gets updated)
    
    [lam0_tf,fval] = fsolve(@cost_minT,lam0_tf_guess,options,t0,x0,xf,a,b,opts_ode);
    
    % Propagate Min-Time Solution
    lam0 = lam0_tf(1:2);
    tf   = lam0_tf(3);
    [t_minT,X_minT] = ode45(@linosc,[t0 tf],[x0'; lam0],opts_ode,a,b);
    
    % Compute Hamiltonian
    H = hamiltonian(t_minT,X_minT,a,b);
    H_minT = H';
    H_minT = H_minT+1;
    
end

%% Analytical Solution

if analytical == true
    
tf      = 12;
A       = [0 1 0 0; -a^2 b 0 -1; 0 0 0 a^2; 0 0 -1 -b];
PHI     = expm(A*tf);   % Phi is the terminal cost we are trying to minimize! 
PHI11   = PHI(1:2,1:2);
PHI12   = PHI(1:2,3:4);
lam0 = inv(PHI12)*(xf' - PHI11*x0');

[t_minU,X_minU] = ode45(@linosc,[t0 tf],[x0'; lam0],opts_ode,a,b);

H = hamiltonian(t_minU,X_minU,a,b);

H_minU = H' + 0.5.*(-X_minU(:,4)).^2;

minU = true;

end
%% Plots
plots(minU,minT,t_minU,X_minU,t_minT,X_minT,H_minU,H_minT);

%% Functions
function Xdot = linosc(t,X,a,b)

x   = X(1:2);
lam = X(3:4);

% State Differential Equations
xdot = [x(2);
        -a^2*x(1) + b*x(2) - lam(2)];
    
% Costate Differential Equations
lam_dot = [a^2*lam(2);
              -(lam(1) + b*lam(2))];
          
Xdot = [xdot; lam_dot];

end

function H = hamiltonian(t,X,a,b)

x   = X(:,1:2);
lam = X(:,3:4);

for i = 1:length(X)
    
    % State Differential Equations
    xdot = [x(i,2);
        -a^2*x(i,1) + b*x(i,2) - lam(i,2)];

    % Hamiltonian
    H(i) = lam(i,:)*xdot;
    
end

end

function err = cost_minU(lam0_guess,t0,tf,x0,xf,a,b,opts_ode)

[t,X] = ode45(@linosc,[t0 tf],[x0'; lam0_guess],opts_ode,a,b);

H = hamiltonian(t,X,a,b);

err = X(end,1:2)' - xf';
err
end

function err = cost_minT(lam0_tf_guess,t0,x0,xf,a,b,opts_ode)

lam0 = lam0_tf_guess(1:2);
tf   = lam0_tf_guess(3);

[t,X] = ode45(@linosc,[t0 tf],[x0'; lam0],opts_ode,a,b);

H = hamiltonian(t,X,a,b);

H = H';

err = [X(end,1:2)'; H(end)+1] - [xf'; 0];

end

function plots(minU,minT,t_minU,X_minU,t_minT,X_minT,H_minU,H_minT)

figure(1)
subplot 211
if minU == true
    plot(t_minU,X_minU(:,1),'b-','Linewidth',2)
end
hold on
if minT == true
    plot(t_minT,X_minT(:,1),'r--','Linewidth',2)
end
if minU == true && minT == true
    legend('Min u','Min Time')
end
ylabel('x_1')
title('States')
subplot 212
if minU == true
    plot(t_minU,X_minU(:,2),'b-','Linewidth',2)
end
hold on
if minT == true
    plot(t_minT,X_minT(:,2),'r--','Linewidth',2)
end
ylabel('x_2')
xlabel('Time')
if minU == true && minT == true
    legend('Fixed Time','Min Time')
end

figure(2)
subplot 211
if minU == true
    plot(t_minU,X_minU(:,3),'b-','Linewidth',2)
end
hold on
if minT == true
    plot(t_minT,X_minT(:,3),'r--','Linewidth',2)
end
if minU == true && minT == true
    legend('Min u','Min Time')
end
ylabel('\lambda_1')
title('Costates')
subplot 212
if minU == true
    plot(t_minU,X_minU(:,4),'b-','Linewidth',2)
end
hold on
if minT == true
    plot(t_minT,X_minT(:,4),'r--','Linewidth',2)
end
ylabel('\lambda_2')
xlabel('Time')
if minU == true && minT == true
    legend('Min u','Min Time')
end

if minU == true
    figure(3)
    subplot 211
    plot(t_minU,H_minU,'b-','Linewidth',2)
    ylabel('H')
    title('Min u')
    subplot 212
    semilogy(t_minU,abs((H_minU-H_minU(1))/H_minU(1)),'b-','Linewidth',2)
    hold on; grid on;
    ylabel('|H-H(t_0)/H(t_0)|')
    xlabel('Time')
end

if minT == true
    figure(4)
    subplot 211
    plot(t_minT,H_minT,'r-','Linewidth',2)
    ylabel('H')
    title('MinT')
    subplot 212
    semilogy(t_minT,abs((H_minT-H_minT(1))),'r-','Linewidth',2)
    hold on; grid on;
    ylabel('|H-H(t_0)/H(t_0)|')
    xlabel('Time')
end


end

