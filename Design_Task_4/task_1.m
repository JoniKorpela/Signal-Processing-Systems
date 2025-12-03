clear;
clc;

L = 3;
M = 4;
f_in = 39e6;
h  = [3, 6, 12, 24, 43, 69, 98, 123, 138, 138, 123, 98, 69, 43, 24, 12, 6, 3];

f1 = 17e6;
f2 = 5e6;

x_n = @(t) sin(2*pi*f1*t) + cos(2*pi*f2*t) + cos(2*pi*0.5e6);


N = 1000;

fs1 = 1;
Ts1 = 1/f_in;
t1 = (1:N) * Ts1;
samples= x_n(t1);
samples_upsampled = zeros(1, N * L);

k=1;
for i=1:N*L
    if (mod(i, L) == 1)
        samples_upsampled(i) = samples(k);
        k = k + 1;
    else
        samples_upsampled(i) = 0;
    end
end


b = h;
a = 1;

filtered_signal = filter(b, a, samples_upsampled);


filtered_signal_downsampled = zeros(1, (N*L) / M);

k=1;
for i=1:N*L
    if(mod(i, M) == 1)
        filtered_signal_downsampled(k) = filtered_signal(i);
        k = k + 1;
    end
end




