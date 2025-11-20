clc;
clear;

%group 2 params for task 2
phi_0 = 8.0;     % Base angle [deg] 
N = 10;          % Number of angles
d_max = 1.0e-3;  % Maximum error

% angle step
delta_phi = 360 / N;

% # of iterations, set arbitrarily
N_iter = 11;

cos_vals = zeros(1,N);
sin_vals = zeros(1,N);
cos_vals_cmp = zeros(1,N);
sin_vals_cmp = zeros(1,N);

for n = 0:N-1
    phi_deg = phi_0 + n*delta_phi;
    phi = deg2rad(phi_deg)

    % phi_deg = -33.75;
    %phi_deg = 100;
    %phi = deg2rad(phi_deg);

    % reduce to +-pi range for correct initialization
    phi = mod(phi + pi, 2*pi) - pi;
    % initialize to +-pi/2 range
    if phi < 0
        dinit = -1;
    else
        dinit = 1;
    end 

    x = 1;
    y = 0;

    x0 = -dinit * y;
    y0 = dinit * x;
    z0 = phi - dinit * pi/2;

    x = x0;
    y = y0;
    z = z0;

    % array for saving rotation dirs
    d_arr = zeros(1, N_iter);

    for k = 0:N_iter-1
    %     if z < 0
    %         if z + atan(2^(-k)) > 0
    %             d = 0;       % don’t rotate
    %         else
    %             d = -1;      % rotate negatively
    %         end
    %     else
    %         if z - atan(2^(-k)) < 0
    %             d = 0;       % don’t rotate
    %         else
    %             d = 1;       % rotate positively
    %         end
    %     end

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

    cos_vals(n+1) = x;
    sin_vals(n+1) = y;
    cos_vals_cmp(n+1) = cosd(phi_deg);
    sin_vals_cmp(n+1) = sind(phi_deg);
    disp('Rotation directions:'); disp(d_arr)
end

disp('cos (CORDIC):'); disp(cos_vals)
disp('cos (true) :'); disp(cos_vals_cmp)
disp('sin (CORDIC):'); disp(sin_vals)
disp('sin (true) :'); disp(sin_vals_cmp)

max_error = max(max(abs(cos_vals_cmp - cos_vals)), max(abs(sin_vals_cmp - sin_vals)))
if max_error >= d_max
    fprintf('Maximum error exceeded!\n');
end

