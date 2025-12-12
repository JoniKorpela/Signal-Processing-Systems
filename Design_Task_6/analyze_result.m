% Script for analyzing result

k = find(bit_error == 1, 1, 'last' )+1;
if k <= 9000
    fprintf('It works!');
    fprintf('\n- no errors after k=%d (%.3f s)',k,k/10000);
    errdiffx = errdiff(k:length(errdiff));
    fprintf('\n - MSE after convergence: %f',mean(errdiffx.^2));
else
    fprintf('It does not work!');
    fprintf('\n- still errors after k = 9000 (0.9 sec)');
    fprintf('\n- adjust parameters (filter length, step size, delay)');
end
fprintf('\n');
