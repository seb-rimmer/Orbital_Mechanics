function circle(x,y,r, varargin)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)


% Check if the optional argument for label is provided
if nargin > 3
    color = varargin{1};
else
    color = '';
end

% Check if the optional argument for color is provided
if nargin > 4
    label = varargin{2};
else
    label = '';
end

ang=0:0.0001:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp, color, 'DisplayName', label);
end

