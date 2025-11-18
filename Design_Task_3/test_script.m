% group 2 coords for task 1
phi = -164.9;    % Angle [deg]
x = 5.25;        % x coordinate
y = 2.75;        % y coordinate

% initialize
phi_rad = deg2rad(phi);

if phi_rad < 0
    dinit = -1;
else
    dinit = 1;
end

x0 = -dinit * y;
y0 = dinit * x;
z0_rad = phi_rad - dinit * pi / 2;

z0 = rad2deg(z0_rad);


fprintf('x0: %.4f, y0: %.4f, z0: %.4f\n', x0, y0, z0_deg);






%group 2 params for task 2
phi_0 = 8.0;     % Base angle [deg] 
N = 10;          % Number of angles
d_max = 1.0e-3;  % Maximum error
