 function Xdot = polar_eom_ocp(t,X,mu)
 
const_and_params;
 
r     = X(1);
theta = X(2);
u     = X(3);
v     = X(4);
lr    = X(5);
lt    = X(6);
lu    = X(7);
lv    = X(8);
 
% Control
U1 = -lu;
U2 = -lv;
 
% State Differential Equations
xdot = [u;
        v/r;
        v^2/r - mu/(r^2) + U1;
        -u*v/r + U2];
    
% Costate Differential Equations
lt = 0;
lam_dot = [lt*(v/r^2) + lu*((v/r)^2 - 2*mu/(r^3)) - lv*(u*v/(r^2));
            0;
            -lr + lv*(v/r);
            -lt/r - lu*(2*v/r) + lv*(u/r)];
 
 
Xdot = [xdot; lam_dot];
 
end

function err = cost_minfuel(lam0_guess,t0,tf,x0,xf,mu,opts_ode)
 
[t,X] = ode45(@polar_eom_ocp,[t0 tf],[x0'; lam0_guess],opts_ode,mu);
 
err = [X(end,1) X(end,6) X(end,3) X(end,4)]' - [xf(1) 0 0 sqrt(mu/xf(1))]'; % final circular orbit
 
end
