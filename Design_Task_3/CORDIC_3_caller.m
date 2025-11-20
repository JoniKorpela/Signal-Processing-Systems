N_iter_max = -1;
N_iter_min = 100;

for x = 0:1:100
    for y = 0:1:100
        N_iter = 1;
        pass = CORDIC_3(x, y, N_iter);
        while pass == false
            N_iter = N_iter + 1;
            if N_iter > N_iter_max
                N_iter_max = N_iter;
            end
            pass = CORDIC_3(x, y, N_iter);
        end
        if N_iter < N_iter_min
            N_iter_min = N_iter;
        end
    end
end

disp(N_iter_max)
disp(N_iter_min)