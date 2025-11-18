clc;
clear;

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


fprintf('x0: %.4f, y0: %.4f, z0: %.4f\n', x0, y0, z0);


%group 2 params for task 2
phi_0 = 8.0;     % Base angle [deg] 
N = 10;          % Number of angles
d_max = 1.0e-3;  % Maximum error

% t = 8 * 36n, n=0, 1,...,9 
%circular system:
% x′(n) = s(n) cos φ
% y′(n) = s(n) sin φ 
% set x = 1 / A_N and y = 0 to evaluate sine/cosine

% angle step
delta_phi = 360 / N;

% # of iterations, set arbitrarily
N_iter = 11;

% CORDIC gain
i = 0:N_iter-1;
A_N = prod(sqrt(1 + 2.^(-2*i)));

cos_vals = zeros(1,N);
sin_vals = zeros(1,N);
cos_vals_cmp = zeros(1,N);
sin_vals_cmp = zeros(1,N);

for n = 0:N-1
    phi_deg = phi_0 + n*delta_phi;
    phi = deg2rad(phi_deg)

    % reduce to +-pi range for correct initialization
    phi = mod(phi + pi, 2*pi) - pi;
    % initialize to +-pi/2 range
    if phi < 0
        dinit = -1;
    else
        dinit = 1;
    end 

    x = 1 / A_N;
    y = 0;

    x0 = -dinit * y;
    y0 = dinit * x;
    z0 = phi - dinit * pi/2

    x = x0;
    y = y0;
    z = z0;

    for k = 0:N_iter-1
        if z < 0
            d = - 1;
        else
            d = 1;
        end

        nx = x - d * y * 2^(-k);
        ny = y + d * x * 2^(-k);
        nz = z - d * atan(2^(-k));
        x = nx;
        y = ny;
        z = nz;
    end

    cos_vals(n+1) = x;
    sin_vals(n+1) = y;
    cos_vals_cmp(n+1) = cosd(phi_deg);
    sin_vals_cmp(n+1) = sind(phi_deg);
end

disp('cos (CORDIC):'); disp(cos_vals)
disp('cos (true) :'); disp(cos_vals_cmp)
disp('sin (CORDIC):'); disp(sin_vals)
disp('sin (true) :'); disp(sin_vals_cmp)

max_error = max(max(abs(cos_vals_cmp - cos_vals)), max(abs(sin_vals_cmp - sin_vals)))
if max_error >= d_max
    fprintf('Maximum error exceeded!\n');
end

[xN, yN, zN, AN] = ucordic(circular, rotation, )
