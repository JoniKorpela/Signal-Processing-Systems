function pass = CORDIC_3(x, y, N_iter)

pass = true;

phi_deg = -33.75;
phi = deg2rad(phi_deg);

% vector to rotate
x_orig = x;
y_orig = y;
z = phi;

% array for saving rotation dirs
d_arr = zeros(1, N_iter);

for k = 0:N_iter-1

    % determine the rotation dir that gets us closest
    c = atan(2^(-k));
    errs = [abs(z - c), abs(z + c), abs(z - 0)];
    [~, idx] = min(errs);
    if idx == 1
        d = 1;
    elseif idx == 2
        d = -1;
    else
        d = 0;
    end
    % store rotation into array for gain computation
    % and display of rotation sequences
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
true_val_x = x_orig * cos(phi) - y_orig * sin(phi);
true_val_y = y_orig * cos(phi) + x_orig * sin(phi);

% compute error
ang_true = atan2(true_val_y, true_val_x);
ang_cordic = atan2(y, x);
ang_error = rad2deg(abs(ang_true - ang_cordic));

% display rotation directions
disp('Rotation directions:'); disp(d_arr)
disp('original angle:'); disp(atan2d(y_orig, x_orig))
disp('true rotated angle:'); disp(rad2deg(ang_true))
disp('CORDIC rotated angle:'); disp(rad2deg(ang_cordic))


error_max = 0.2; %deg
disp('Angular error:'); fprintf('%.8f\n', ang_error);
if max(ang_error) > error_max
    fprintf('Maximum error exceeded!\n\n');
    pass = false;
end

