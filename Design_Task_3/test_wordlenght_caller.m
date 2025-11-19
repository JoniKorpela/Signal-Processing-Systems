
clc;
clear;

% sweep
xy_p_values = 10:20;   % total bits for x,y
xy_n_values = 8:16;    % fractional bits for x,y

z_p_values  = 10:20;   % total bits for angle z
z_n_values  = 8:16;    % fractional bits for angle z

best = [];  % store best configuration

for p1 = xy_p_values
    for n1 = xy_n_values

        % skip illegal combinations
        if n1 >= p1, continue; end

        for p2 = z_p_values
            for n2 = z_n_values

                if n2 >= p2, continue; end

                % build struct
                fixptpar.xy_p = p1;
                fixptpar.xy_n = n1;
                fixptpar.z_p  = p2;
                fixptpar.z_n  = n2;

                fprintf('\nTesting xy=%d.%d  z=%d.%d\n', p1,n1,p2,n2);

                % Run the test
                    if test_wordlength(fixptpar);
                        fprintf('Passed: xy=%d.%d  z=%d.%d\n', p1,n1,p2,n2);
                        best = [best; p1 n1 p2 n2];
                    end;
            end
        end
    end
end

disp('===============================');
disp('All passing configurations:');
disp('xy_p   xy_n   z_p   z_n');
disp(best);
disp('===============================');