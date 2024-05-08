r1 = 1;
r2 = 5.2;
theta = 147;

circle(0, 0, r1)
grid on 
hold on
circle(0, 0, r2)
abs_lim = r2*1.1;
xlim([-abs_lim abs_lim])
ylim([-abs_lim abs_lim])

% first focii at origin, vacant focci tbd
f1 = [0 0]';
f2 = [0 0]';

p1 = [r1 0];
theta = theta * pi/180;
p2 = [r2*cos(theta) r2*sin(theta)];

% size of vectors for radii and space triangle
% r1 = norm(p1);
% r2 = norm(p2);
p_12 = p2-p1;
r_12 = norm(p_12);

% Min energy transfer ellipse parameters
s = 0.5 * (r1 + r2 + r_12);
a_min = s/2;

% calculating position of vacant focii
p1_fvac = p_12 * (2*a_min - r1)/norm(p_12);
p1_fvac_norm = norm(p1_fvac);
f_fvac = p1_fvac + [1 0];
p1_fvac_origin = p1_fvac + [1 0];

e = norm(f_fvac)/(2*a_min);
tilt_ellipse = atan2(p1_fvac_origin(2), p1_fvac_origin(1));
b = a_min * sqrt(1 - e^2);

% Transfer ellipse points
actual_center = p1_fvac_origin/2;
angle = tilt_ellipse;
theta = linspace(0, 2*pi, 10000); % angle from 0 to 2*pi
a = a_min;
x_y = [a * cos(theta); b * sin(theta);];
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
x_y = R * x_y;
x_y =  x_y + actual_center';

% r1
plot([0 p1(1)], [0 p1(2)], 'k')
% r2
plot([0 p2(1)], [0 p2(2)], 'k')
% Vertice of space triangle
plot([1 1+p_12(1)], [0 p_12(2)]);

% Line from r1 to vacant focii
plot([1 p1_fvac_origin(1)], [0 p1_fvac_origin(2)], 'LineWidth', 2);

% Plot transfer ellipse
plot(x_y(1,:), x_y(2,:));
scatter(actual_center(1), actual_center(2))

% Plot focii
scatter(0, 0, '*', 'r') % F
scatter(p1_fvac_origin(1), p1_fvac_origin(2), '*', 'r') % F*

axis equal