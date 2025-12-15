
function main

    % Set initial L to 100 instead of 400, since clearly 400 is too much 
    % for a system with taps < 100
    L = 100;
    % set step size according to 0 < μ < 1/(LPu). Start from (1/(LPu)) / 2
    % grab u from T2data for P_u calculation
    [u, ~] = T2data(2);
    P_u = mean(u.^2);
    step_size = (1/(L*P_u))/2;
    
    % coarse search -> fine search

    % Coarse search for delay D
    [D_low, D_high, D_coarse] = find_D(10, 10, 200);
    % Fine search around the coarse value of D
    [~, ~, D] = find_D(D_coarse - 10, 1, D_coarse + 10);

    % Coarse search for filter length L
    [L_low, L_high, L_coarse] = find_L(10, 10, 300);
    % Fine search around the coarse value of L
    [~, ~, L] = find_L(max(1, L_coarse - 10), 1, L_coarse + 10);

    % Coarse search for step size
    [sz_low, sz_high, step_size_coarse, ~, ~] = find_step_size((1/(L*P_u))/10, (1/(L*P_u))/100, 1/(L*P_u));
    delta = 0.2 * step_size_coarse; % Fine-tuning range for step size
    start_val = max(step_size_coarse - delta, 0);
    end_val = min(step_size_coarse + delta, 1/(L*P_u)); % mu max
    step = delta / 10;
    % Fine search around coarse value of step size
    [~, ~, step_size, cvrg_best, mse] = find_step_size(start_val, step, end_val);
    
    % Print best parameters
    fprintf("\nD = %d\n", D);
    fprintf("L = %d\n", L);
    fprintf("step_size = %.8f\n", step_size);
    fprintf("Converges at k = %d (%.3f s)\n", cvrg_best, cvrg_best/10000);
    fprintf("MSE after convergence: %.6f\n\n", mse);

    % Print ranges of passing values
    fprintf("Range of passing D's: [%d, %d]\n", D_low, D_high);
    fprintf("Range of passing L's: [%d, %d]\n", L_low, L_high);
    fprintf("Range of passing step sizes: [%.6f, %.6f]\n", sz_low, sz_high);

    % Run exsim with best parameters
    exsim_run(L, D, step_size);


    % find smallest D that converges within constraint. Also return lowest
    % and highest passable D
    function [d_low, d_high, d_out] = find_D(start_val, step, end_val)
        d_low = Inf;
        d_high = -Inf;
        d_out = NaN;
        found = false;
        for d = start_val:step:end_val
            disp(d)
            [~, pass, ~, ~] = exsim_run(L, d, step_size);
            if pass == true
                if found == false
                    d_out = d; % lowest passing D is found, don't overwrite
                    found = true; 
                end
                if d > d_high
                    d_high = d;
                end
                if d < d_low
                    d_low = d;
                end
            end
        end
    end
    % select the L that converges fastest. Return lowest and highest
    % passable L
    function [l_low, l_high, l_out] = find_L(start_val, step, end_val)
        best = Inf;
        l_out = NaN;
        l_low = Inf;
        l_high = -Inf;
        for l = start_val:step:end_val
            disp(l)
            step_size = (1/(l*P_u))/2;
            [~, pass, cvrg, ~] = exsim_run(l, D, step_size);
            if pass == true
                if cvrg < best
                    best = cvrg;
                    l_out = l;
                end
                if l > l_high
                    l_high = l; 
                end
                if l < l_low
                    l_low = l; 
                end
            end
        end
    end
    % Select the sz that converges fastest. Return convergence step and mse
    % also. Return lowest and highest passing values.
    function [lowest_sz, highest_sz, sz_out, cvrg_out, mse_out] = find_step_size(start_val, step, end_val)
        best = Inf;
        sz_out = NaN;
        lowest_sz = Inf;
        highest_sz = -Inf;
        for sz = start_val:step:end_val
            disp(sz)
            [~, pass, cvrg, MSE] = exsim_run(L, D, sz);
            if pass == true
                if cvrg < best
                    best = cvrg;
                    cvrg_out = cvrg;
                    sz_out = sz;
                    mse_out = MSE;
                end
                if sz < lowest_sz
                    lowest_sz = sz;
                end
                if sz > highest_sz
                    highest_sz = sz;
                end
            end
        end
    end
end

main()