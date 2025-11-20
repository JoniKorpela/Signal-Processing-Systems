clc;
clear;

%group 2 params for task 2
%CORDIC-3
N = 10;
% # of iterations, set arbitrarily
N_iter = 2;

vals_x = zeros(1,N);
vals_y = zeros(1,N);
true_vals_x = zeros(1,N);
true_vals_y = zeros(1,N);
ang_errors = zeros(1,N);

for n = 0:N-1
    % phi_deg = phi_0 + n*delta_phi;


    phi_deg = -33.75;
    phi = deg2rad(phi_deg);

    % reduce to +-pi range for correct initialization
    phi = mod(phi + pi, 2*pi) - pi;
    % initialize to +-pi/2 range
    if phi < 0
        dinit = -1;
    else
        dinit = 1;
    end 

    % vector to rotate
    x_orig = 10;
    y_orig = 1;

    x0 = -dinit * y_orig;
    y0 = dinit * x_orig;
    z0 = phi - dinit * pi/2;

    x = x0;
    y = y0;
    z = z0;

    % array for saving rotation dirs
    d_arr = zeros(1, N_iter);

    for k = 0:N_iter-1

        % determine the rotation dir that gets us closest
        angles = [atan(2^(-k)), -atan(2^(-k)), 0];
        [~, idx] = min( abs( z - angles ) );
        if idx == 1
            d = 1;
        elseif idx == 2
            d = -1;
        else
            d = 0;
        end

        % store rotation into aray for gain computation
        d_arr(k + 1) = d;

        nx = x - d * y * 2^(-k);
        ny = y + d * x * 2^(-k);
        nz = z - d * atan(2^(-k));
        x = nx;
        y = ny;
        z = nz;
    end

    % CORDIC gain
    A_N = 1;
    for k = 0:N_iter-1
        if d_arr(k+1) ~= 0
            A_N = A_N * sqrt(1 + 2^(-2*k));
        end
    end
    x = x / A_N ;
    y = y / A_N ;

    % store true and computed vector values
    vals_x(n+1) = x;
    vals_y(n+1) = y;
    true_vals_x(n+1) = x_orig * cos(phi) - y_orig * sin(phi);
    true_vals_y(n+1) = y_orig * cos(phi) + x_orig * sin(phi);

    % compute error
    ang_true = atan2(true_vals_y(n+1), true_vals_x(n+1));
    ang_cordic = atan2(y, x);
    ang_errors(n+1) = abs(ang_true - ang_cordic);

    % display rotation directions
    disp('Rotation directions:'); disp(d_arr)

end
disp('x (CORDIC):'); disp(vals_x)
disp('x (true) :'); disp(true_vals_x)
disp('y (CORDIC):'); disp(vals_y)
disp('y (true) :'); disp(true_vals_y)

error_max = 0.2; %deg
disp('Maximum angular error:'); disp(max(ang_errors));
if max(ang_errors) >= error_max
    fprintf('Maximum error exceeded!\n');
end

