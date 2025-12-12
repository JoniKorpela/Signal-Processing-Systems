clc;
clear;

h = [-2, 9, -24, 33, -24, 8, -1];
u = [1, 1, 0, 2, -2, 0, 1, -1, 1, 2, -2, 1, 0, 1,...
    -1, 1, 2, -1, 1, 2, 0, 1, -1, 2, -1, 0, -2, -1, 2, -1, -2, -1];

% Doing the linear convolution just to compare results
y_linear = conv(h, u);

N_h = 7;        % impulse response length
N_u = 6;        % section length
N = 16;         % FFT length. N_u + N_h - 1 -->> round to power of 2

h_pad = [h, zeros(1, N - N_h)];

% convolution output, initialize empty and concatenate
y_ols = [];

u = [zeros(1,N_h-1), u];  % prepend zeros
% compute number of blocks
num_sections = ceil((length(u) - (N_h-1)) / (N - N_h + 1));
total_length = num_sections*(N - N_h + 1) + N_h - 1; 
% append zeroes to make final block full length
u = [u, zeros(1, total_length - length(u))]; 

for k = 0:num_sections-1
    start_idx = k*(N - N_h + 1) + 1;
    u_seg = u(start_idx : start_idx + N - 1);
    y_seg = cyclic_conv(u_seg, h_pad, N);
    y_valid = y_seg(N_h:end);      % discard first Nh samples
    y_ols = [y_ols, y_valid];      % concatenate
end

disp('OLS convolved data:')
disp(y_ols)

figure(1);
hold on;
stem(y_ols, 'filled', Color='g');
stem(y_linear, Color='r');
legend('OLS convolution', 'conv()')
hold off;
