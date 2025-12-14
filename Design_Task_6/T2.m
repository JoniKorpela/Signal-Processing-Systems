
function main
% set initial L to 100 instead of 400, since clearly 400 is too much 
% for a system with taps < 100
L = 100; 
D = 40;
% set step size according to 0 < μ < 1/(LPu). Start from (1/(LPu)) / 2
% grab u from T2data
[u, ~] = T2data(2);
P_u = mean(u.^2);
step_size = (1/(L*P_u))/2;

% coarse search -> fine search
D_coarse = find_D(10, 10, 200);
disp(D_coarse)
D = find_D(D_coarse - 10, 1, D_coarse + 10);
disp(D)
L_coarse = find_L(10, 10, 200);
L = find_L(max(1, L_coarse - 10), 1, L_coarse + 10);
step_size_coarse = find_step_size((1/(L*P_u))/10, (1/(L*P_u))/100, 1/(L*P_u));
delta = 0.2 * step_size_coarse;
start_val = max(step_size_coarse - delta, 0);
end_val = min(step_size_coarse + delta, 1/(L*P_u)); % mu max
step = delta / 10;
[step_size, cvrg_best, mse] = find_step_size(start_val, step, end_val);

fprintf("D = %d\n", D);
fprintf("L = %d\n", L);
fprintf("step_size = %.8f\n", step_size);
fprintf("Converges at k = %d (%.3f s)", cvrg_best, cvrg_best/10000);
fprintf("MSE after convergence: %.6f", mse);

    % find smallest D that converges within constraint
    function [d_out] = find_D(start_val, step, end_val)
        d_out = NaN;
        for d = start_val:step:end_val
            disp(d)
            [~, pass, ~, ~] = exsim_run(L, d, step_size);
            if pass == true
                d_out = d;
                break;
            end
        end
    end
    % select the L that converges fastest
    function [l_out] = find_L(start_val, step, end_val)
        best = Inf;
        l_out = NaN;
        for l = start_val:step:end_val
            disp(l)
            step_size = (1/(l*P_u))/2;
            [~, pass, cvrg, ~] = exsim_run(l, D, step_size);
            if pass == true
                if cvrg < best
                    best = cvrg;
                    l_out = l;
                end
            end
        end
    end
% select the sz that converges fastest. return convergence step and mse
% also
    function [sz_out, cvrg_out, mse_out] = find_step_size(start_val, step, end_val)
        best = Inf;
        sz_out = NaN;
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
            end
        end
    end
end


main()