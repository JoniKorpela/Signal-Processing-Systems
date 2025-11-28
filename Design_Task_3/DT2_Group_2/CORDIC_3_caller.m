clc;
clear;

N_iter_max = 0;
N_iter_min = 100;

for x = 0:0.1:10
    for y = 0:0.1:10
        N_iter = 5;
        pass = CORDIC_3(x, y, N_iter);
        while pass == false
            N_iter = N_iter + 1;
            pass = CORDIC_3(x, y, N_iter);
        end
        if N_iter < N_iter_min
            N_iter_min = N_iter;
        end
        if N_iter > N_iter_max
            N_iter_max = N_iter;
        end
    end
end

disp('Lowest N_iter to pass: '); disp(N_iter_min);
disp('Minimum N_iter required to pass all'); disp(N_iter_max);


