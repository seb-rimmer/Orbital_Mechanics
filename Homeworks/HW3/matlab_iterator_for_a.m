% Lambert Solver for HW3 Problem_1 - Earth-Apophis Intercept

clear all; 
close all; clc;

a = fzero(@lambert, 1.1)

function f = lambert(a)

s  = 1.8381;        % Pre-calculated space-triangle semi-perimeter
c  = 1.6007;        % Pre-calculated chord length
tf = (410/365.25) * 2 * pi;     % ToF of 410 days

% Updated value for alpha because our t_f is greater than t_m 
alpha = 2*pi - 2*asin(sqrt(s/(2*a)));
beta  = - 2*asin(sqrt((s-c)/(2*a)));

f = tf-(a^(3/2))*(alpha-beta-sin(alpha)+sin(beta));

end
