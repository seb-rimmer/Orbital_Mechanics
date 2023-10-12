% Lambert Solver for HW3 Problem_1 - Earth-Apophis Intercept

clear all;
close all; clc;

function f = lambert(a)

s  = 6.13;        % Pre-calculated space-triangle semi-perimeter
c  = 6.06;        % Pre-calculated chord length
tf = (524/365.25) * 2 * pi;     % ToF of 410 days

% Updated value for alpha because our t_f is greater than t_m
alpha = 2*asin(sqrt(s/(2*a)));
beta  = 2*asin(sqrt((s-c)/(2*a)));

f = tf-(a^(3/2))*(alpha-beta-sin(alpha)+sin(beta));

end

a = fzero(@lambert, 5)