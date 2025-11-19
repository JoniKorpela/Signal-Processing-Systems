
%group 2 params for task 2
phi_0 = 8.0;     % Base angle [deg] 
N = 10;          % Number of angles
d_max = 1.0e-3;  % Maximum error

% angle step
delta_phi = 360 / N;

% # of iterations to try
N_iter_max = 20;

% fixed point params
% p1 = 20;   % xy_p
% n1 = 14;  % xy_n
% p2 = 20;   % z_p
% n2 = 14;  % z_n
% fixptpar.xy_p = p1;
% fixptpar.xy_n = n1;
% fixptpar.z_p  = p2;
% fixptpar.z_n  = n2;
gaincomp = true;
override = false;
verbose = false;

for n = 0:N-1
    
    phi_deg = phi_0 + n*delta_phi;
    phi = deg2rad(phi_deg);
    phi = mod(phi + pi, 2*pi) - pi;
    
    
    if phi < 0
        d_init = -1;
    else
        d_init = 1;
    end
    x0 = 0;
    y0 = d_init;
    z0 = phi - d_init * pi / 2;
    
    % true values
    cos_true(n+1) = cos(phi);
    sin_true(n+1) = sin(phi);
    
    found_ok = false;
    for k = 1:N_iter_max
        [xN, yN, ~, ~] = ucordic('circular','rotation', x0, y0, z0, k, gaincomp, override, fixptpar, verbose);
        cos_cordic(n+1) = xN;
        sin_cordic(n+1) = yN;
    
        max_err = max( abs(cos_true(n+1)-cos_cordic(n+1)), abs(sin_true(n+1)-sin_cordic(n+1)) );
        if max_err <= d_max
            fprintf('Angle %.1f°: reached error %.3g with k=%d iterations\n', phi_deg, max_err, k);
            found_ok = true;
            break;
        end
    end
    
    if ~found_ok
        fprintf('Angle %.1f°: did NOT reach error %.3g with up to %d iterations (best error = %.3g)\n', ...
                phi_deg, d_max, N_iter_max, max_err);
    end
end

disp('cos (CORDIC):'); disp(cos_cordic);
disp('cos (true) :'); disp(cos_true);
disp('sin (CORDIC):'); disp(sin_cordic);
disp('sin (true) :'); disp(sin_true);

