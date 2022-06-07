% Robyn Woollands
% AE 402 - Lec 11 Example

clear all; close all; clc;


a = fzero(@lambert,5)

function f = lambert(a)

s  = 6.13;
c  = 6.06;
tf = (524/365.25)*2*pi;

alpha = 2*asin(sqrt(s/(2*a)));
beta  = 2*asin(sqrt((s-c)/(2*a)));

f     = tf-(a^(3/2))*(alpha-beta-sin(alpha)+sin(beta));

end
