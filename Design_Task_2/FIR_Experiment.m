% FIR_EXPERIMENT Function for performing experiments in the exercise H2
%
% << Syntax >>
%    [y_fixpt,y_float,y_error] = FIR_Experiment(input_mode,a,bitspec,floatmode)
% 
% << Arguments >>
%    input_mode    Select the input signal : 'sine', 'rand' or 'intro' 
%                   - the first two are used in the design task 3 !!
%    a             Unquantized FIR filter coefficients
%    bitspec       Word/fraction length specification
%    y_fixpt       Fixed-point filtering result
%    y_float       Floating-point filtering result
%    y_error       Their difference 
%
% << Examples >>
%    % Running the example provided in intro2.pdf, page 18
%    Hd = x_filter; a = Hd.Numerator; bitspec = [8,8,7,0,6,4];
%    [y_fixpt,y_float,y_error] = FIR_Experiment('intro',a,bitspec);

function [y_fixpt,y_float,y_error] = FIR_Experiment(input_mode,a,bitspec)

% Using FFT of length N to analyze simulation results
N = 256;

% The number of extra data points: only the last N points
% of the filter output will be used to display results
extra = length(a);
datalen = N + extra;

% Generate the input data
 if strcmp(input_mode,'sine')

    % This signal is used as the first choice in Task 3(b)
    % The sampling frequency is 30 kHz. The sum of three sinusoids 
    % (3.5 kHz, 7.5 kHz, and 8.5 kHz) is sampled.
    
    % Sum of three sine waves
    % -----------------------
    Fs = 30;

    % Three sine waves
    frequency1 = 3.5/Fs; 
    frequency2 = 7.5/Fs;
    frequency3 = 8.5/Fs;
    t = (0:(datalen-1));
    x = 0.35 * sin(2 * pi * frequency1 * t) ...
      + 0.3 * sin(2 * pi * frequency2 * t) ...
      + 0.35 * sin(2 * pi * frequency3 * t);
 
elseif strcmp(input_mode,'rand')

    % This signal is used as the second choice in Task 3(b).
    % Random data, uniform distribution over [-1,1], no correlation
    % The signal is considered to be sampled at 25 kHz.
  
    s = RandStream('mt19937ar','Seed',1275);
    x = 2 * rand(s,1,datalen) - 1;

elseif strcmp(input_mode,'intro')

    % This signal is used in the example of the introduction document
    % The sampling frequency is 40 kHz. The sum of two sinusoids (8 kHz, 14
    % kHz) is sampled.

    Fs = 40;

    % Three sine waves
    frequency1 = 8/Fs; 
    frequency2 = 14/Fs;
    t = (0:(datalen-1));
    x = 0.5 * (sin(2 * pi * frequency1 * t) ...
               + sin(2 * pi * frequency2 * t));

else
    
    error('invalid input_mode parameter');
    
end

if length(bitspec) ~= 6
    error('6 bitspec parameters expected');
end

p = bitspec(1);
p_c = bitspec(2);
n_c = bitspec(3);
r = bitspec(4);
g = bitspec(5);
n = bitspec(6);

fprintf(1,'I/O word length: %d\n',p);
fprintf(1,'MAC input format: s%d.%d\n',p,p-1);
fprintf(1,'Coefficient format: s%d.%d\n',p_c,n_c);
fprintf(1,'Product pipeline register format: s%d.%d\n',p+p_c-r,p-1+n_c-r);
fprintf(1,'Accumulator register format: s%d.%d\n',p+p_c-r+g,p-1+n_c-r);
fprintf(1,'MAC output format: s%d.%d\n',p,n);

% With this scaling some overflow warnings are avoided
scale = (2^p)/(2^p-1);
x = x / scale;

% Setup input arguments for FIR_MAC_Sim
qxp = [p,p-1,1];
qap = [p_c,n_c,1];
qpp = [p+p_c-r,p-1+n_c-r];
gbits = g;
qop = [p,n];

% Perform fixed-point and floating-point simulations
h = FIR_MAC_Simulator(qxp,qap,qpp,gbits,qop);
float_s = h.simulate(1,x,a); % Double precision result
fixpt_s = h.simulate(0,x,a); % Fixed-point simulation result
diff_s = h.computeDiff(float_s);

% Visualize the result in time-domain
indices =  (datalen - N + 1) : datalen;
x = x(indices);
y_fixpt = fixpt_s.out(indices);
y_float = float_s.out(indices);
y_error = diff_s.out(indices);

figure(1);

t = 1:N;
subplot(411); 
plot(t,x); 
set(gca,'xlim',[0 N+1]);
title('Input signal','fontsize',14); set(gca,'fontsize',14);

subplot(412); 
plot(t,y_float); 
set(gca,'xlim',[0 N+1]);
title('Floating-point filtering result','fontsize',14); set(gca,'fontsize',14);

subplot(413); 
plot(t,y_fixpt); 
set(gca,'xlim',[0 N+1]);
title('Fixed-point filtering result','fontsize',14); set(gca,'fontsize',14);

subplot(414); 
plot(t,y_error);
set(gca,'xlim',[0 N+1]);
title('Difference of results','fontsize',14); set(gca,'fontsize',14);
