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

% pad h and u to length N_u + N_h + N_f - 1 = 16

L_valid = N_u + N_h - 1; % valid samples per block = 12 
num_sections = ceil(length(u)/N_u);

% pad tail of u to complete last block
u = [u, zeros(1, num_sections * N_u - length(u))];

% convolution output. Length = L_u + L_h - 1
y_ola = zeros(1, length(u) + length(h) - 1);

% pad h to N bits = 16
h_pad = [h, zeros(1, N - N_h)];

for k = 0:num_sections - 1
    % take section
    u_seg = u(k*N_u + 1 : k*N_u + N_u);
    
    % pad to 16
    u_pad = [u_seg, zeros(1, N - N_u)];

    % cyclic convolution
    y_seg = cyclic_conv(u_pad, h_pad, N);

    % keep the valid 12 samples
    y_valid = y_seg(1:L_valid);

    % indexes for valid OLA segment
    start_idx = k * N_u + 1;
    end_idx = start_idx + L_valid - 1;
   
    if end_idx > length(y_ola)
        end_idx = length(y_ola);
        y_valid = y_valid(1 : end_idx - start_idx + 1);
    end

    % write to output array
    y_ola(start_idx : end_idx) = y_ola(start_idx : end_idx) + y_valid;
end

disp('OLA convolved data:')
disp(y_ola)

figure(1);
hold on;
stem(y_ola, 'filled', Color='g');
stem(y_linear, Color='r');
legend('OLA convolution', 'conv()')
hold off;


