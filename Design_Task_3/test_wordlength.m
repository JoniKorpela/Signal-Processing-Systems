function pass = test_wordlength(fixptpar)

% group 2 params
phi_0 = 8.0;
N = 10;
d_max = 1.0e-3;

delta_phi = 360 / N;
N_iter = 11; % determined through empirical testing

gaincomp = true;
override = false;
verbose = false;

pass = true;   % success unless violation

for n = 0:N-1

    phi_deg = phi_0 + n * delta_phi;
    phi = deg2rad(phi_deg);
    phi = mod(phi + pi, 2*pi) - pi;

    if phi < 0
        d_init = -1;
    else
        d_init = 1;
    end

    x0 = 0;
    y0 = d_init;
    z0 = phi - d_init * pi/2;

    cos_true = cos(phi);
    sin_true = sin(phi);

    try
        [xN, yN, ~, ~] = ucordic('circular','rotation', ...
                                 x0, y0, z0, N_iter, ...
                                 gaincomp, override, fixptpar, verbose);
    catch
        % error -> fail
        pass = false;
        return;
    end

    max_err = max(abs(cos_true - xN), abs(sin_true - yN));

    if max_err > d_max % max error exceeded -> fail
        pass = false;
        return;
    end

end

% if loop completes, everything passed
end