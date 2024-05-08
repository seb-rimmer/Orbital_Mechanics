%% lamberts_visual_choose_a

% Purpose of this function is to visualise the possible trajectories
% between two points at r1 and r2 separated by theta using elliptical
% transfer. Lamberts problem is not solved, time of slight not considered,
% purely a visual tool for possible transfers.

r1 = 1;
r2 = 5.2;
theta = 147;

% chosen value for a
a_chosen = 4;

circle(0, 0, r1, 'r', 'R1')
grid on 
hold on
circle(0, 0, r2, 'b', 'R2')
abs_lim = r2*2;
xlim([-abs_lim abs_lim])
ylim([-abs_lim abs_lim])

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

%% Focii for two ellipses for a_chosen
syms x y
eq1 = (x - p1(1))^2 + (y - p1(2))^2 == (2*a_chosen - r1)^2;
eq2 = (x - p2(1))^2 + (y - p2(2))^2 == (2*a_chosen - r2)^2;
sol = solve([eq1, eq2], [x, y]);

f1 = [double(sol.x(1)) double(sol.y(1)) ];
f2 = [double(sol.x(2)) double(sol.y(2)) ];

e1 = norm(f1)/(2*a_chosen);
tilt_ellipse1 = atan2(f1(2), f1(1));
b1 = a_chosen * sqrt(1 - e1^2);

e2 = norm(f2)/(2*a_chosen);
tilt_ellipse2 = atan2(f2(2), f2(1));
b2 = a_chosen * sqrt(1 - e2^2);

% Min energy transfer ellipse points
actual_center1 = f1/2;
angle1 = tilt_ellipse1;
theta1 = linspace(0, 2*pi, 10000); % angle from 0 to 2*pi
x_y1 = [a_chosen * cos(theta1); b1 * sin(theta1);];
R = [cos(angle1) -sin(angle1); sin(angle1) cos(angle1)];
x_y1 = R * x_y1;
x_y1 =  x_y1 + actual_center1';

% Min energy transfer ellipse points
actual_center2 = f2/2;
angle2 = tilt_ellipse2;
theta2 = linspace(0, 2*pi, 10000); % angle from 0 to 2*pi
x_y2 = [a_chosen * cos(theta2); b2 * sin(theta2);];
R = [cos(angle2) -sin(angle2); sin(angle2) cos(angle2)];
x_y2 = R * x_y2;
x_y2 =  x_y2 + actual_center2';

%% Calculating position of vacant focii for a_min ellipse
p1_fvac = p_12 * (2*a_min - r1)/norm(p_12);
p1_fvac_norm = norm(p1_fvac);
f_fvac = p1_fvac + [1 0];
p1_fvac_origin = p1_fvac + [1 0];

e = norm(f_fvac)/(2*a_min);
tilt_ellipse = atan2(p1_fvac_origin(2), p1_fvac_origin(1));
b = a_min * sqrt(1 - e^2);

% Min energy transfer ellipse points
actual_center = p1_fvac_origin/2;
angle = tilt_ellipse;
theta = linspace(0, 2*pi, 10000); % angle from 0 to 2*pi
a = a_min;
x_y = [a * cos(theta); b * sin(theta);];
R = [cos(angle) -sin(angle); sin(angle) cos(angle)];
x_y = R * x_y;
x_y =  x_y + actual_center';

%% Plotting

% r1
plot([0 p1(1)], [0 p1(2)], 'k')
% r2
plot([0 p2(1)], [0 p2(2)], 'k')
% Vertice of space triangle
% plot([1 1+p_12(1)], [0 p_12(2)]);

% Line from r1 to vacant focii
% plot([1 p1_fvac_origin(1)], [0 p1_fvac_origin(2)], 'LineWidth', 2);


% Plot circle centered at P1
circle(p1(1), p1(2), 2*a_chosen - r1, 'k--', 'Circle for p1');

% Plot circle centered at P2
circle(p2(1), p2(2), 2*a_chosen - r2, 'k--', 'Circle for p2');

% Plot minimum energy transfer ellipse
plot(x_y(1,:), x_y(2,:), 'magenta', 'DisplayName', 'Minimum energy transfer ellipse');

% Plotting ellipses for trajectories
dark_green = [0, 0.5, 0]; % RGB values for dark green
plot(x_y1(1,:), x_y1(2,:), 'Color', dark_green, 'LineStyle', '--', 'DisplayName', 'Transfer ellipse 1 - a_chosen');
plot(x_y2(1,:), x_y2(2,:), 'Color', dark_green, 'LineStyle', '--', 'DisplayName', 'Transfer ellipse 2 - a_chosen');

% scatter(actual_center(1), actual_center(2))

% Plot focii
% scatter(0, 0, '*', 'r') % F
% scatter(p1_fvac_origin(1), p1_fvac_origin(2), '*', 'r') % F*

legend('show');
axis equal