% b 
clear;

% Number of stages; Downsampling factor; Diff delay factor
N = 4;              R = 3;               M = 2;
p = 9;  % wordlen
n = 7;  % fraction bits
x = 20; % output wordlen


% Uncomment to simulate with constant 1

% sysfun = cic('decim',N,R,M);
% arr = ones(1, 1000)
% [output,iplot] = sysfun(arr,0,p,x,7);
% figure(1); iplot(1)
% figure(2); iplot(2)
% gain = (R * M)^N;
% figure(3); plot(output/gain)


% Uncomment to simulate with sine

sysfun = cic('decim',N,R,M);
idx = 0:0.1:99.9
arr = sin(idx);
[output,iplot] = sysfun(arr,0,p,x,7);
figure(1); iplot(1)
figure(2); iplot(2)
gain = (R * M)^N;
figure(3); plot(output/gain)