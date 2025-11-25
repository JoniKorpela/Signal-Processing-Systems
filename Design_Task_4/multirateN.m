% MULTIRATEN  Tool for N-stage multirate design experiment
%
% << Syntax >>
%    multirateN(H0spec,Ms,HK_Ap,HK_As,HA_fs,HA_Ap,HA_As)
%
% << Arguments >>
%    H0spec    Specification of the target FIR : structure with fields
%              .F0    Sampling frequency           [kHz]
%              .fp    Passband cutoff frequency    [kHz]
%              .fs    Stopband cutoff frequency    [kHz]
%              .Ap    Maximum passband ripple      [dB]
%              .As    Minimum stopband attenuation [dB] (positive value)
%    M_v       Down/upsampling factors  : vector of N elements
%    HK_Ap     Maximum passband ripple of the kernel filter [dB]
%    HK_As     Minimum stopband attenuation of the kernel [dB]
%    HA_fs_v   Stopband cutoff frequency of antialias/anti-image filters : N-element vector
%    HA_Ap_v   Passband ripple of antialias/anti-image filters : N-element vector
%    HA_As_v   Minimum stopband attenuation of antialias/anti-image filters : N-element vector
%              - N >= 1 is the number of decimation/interpolation stages
%
% << Description >>
%    The tool can be used to experiment multirate lowpass FIR design, where
%    one or more decimation/interpolation stages are used. Specification of
%    the filters is provided as an input, and the "fvtool" is then used to verify
%    the multirate design. The function uses "firpm" for the filter design.
%
% << Examples >>
%    % --------------------------------------
%    % Example 2 from the background material
%    % --------------------------------------
%    H0spec.F0 = 48; H0spec.fp = 1.0; H0spec.fs = 1.2; 
%    H0spec.Ap = 0.1; H0spec.As = 60;
%    M_v = [4];                   % Single stage, M = 4
%    HK_Ap = 0.03;                % Kernel filter
%    HK_As = 61;
%    HA_fs_v = [10.8];            % Stopband cutoff of decimation filter
%    HA_Ap_v = [0.02];            % Passband ripples
%    HA_As_v = [65];              % Stopband attenuations
%    multirateN(H0spec,M_v,HK_Ap,HK_As,HA_fs_v,HA_Ap_v,HA_As_v)
%    % --------------------------------------
%    % Example 3 from the background material
%    % --------------------------------------
%    H0spec.F0 = 48; H0spec.fp = 1.0; H0spec.fs = 1.2; 
%    H0spec.Ap = 0.1; H0spec.As = 60;
%    M_v = [2,2];                   % For both stages the down/upsampling factor is 2 
%    HK_Ap = 0.03;                  
%    HK_As = 61;
%    HA_fs_v = [22.8,10.8];         
%    HA_Ap_v = [0.007,0.007];       
%    HA_As_v = [78,70];             
%    multirateN(H0spec,M_v,HK_Ap,HK_As,HA_fs_v,HA_Ap_v,HA_As_v)

% $Id: multirateN.m,v 1.1 2013/11/01 06:43:40 psangi Exp $

function multirateN(H0spec,M_v,HK_Ap,HK_As,HA_fs_v,HA_Ap_v,HA_As_v)
   
    Nstages = length(M_v);
    if length(HA_fs_v) ~= Nstages
        error('length of 4th arg (HA_fs_v) must be equal to the length of 2nd (M_v).\n');
    end
    if length(HA_Ap_v) ~= Nstages
        error('length of 5th arg (HA_Ap_v) must be equal to the length of 2nd (M_v).\n');
    end
    if length(HA_As_v) ~= Nstages
        error('length of 6th arg (HA_As_v) must be equal to the length of 2nd (M_v).\n');
    end
    
    % Single-rate filter
    Hd0 = design_filter(H0spec.F0*1000,H0spec.fp*1000,H0spec.fs*1000,H0spec.Ap,H0spec.As);
        
    % Multirate kernel filter
    M = prod(M_v);
    FK = H0spec.F0 / M; % Sampling frequency of the kernel stage
    HdK = design_filter(FK*1000,H0spec.fp*1000,H0spec.fs*1000,HK_Ap,HK_As);
    
    % Multirate decimation / interpolation filters
    HdA_v = cell(Nstages);
    Min = 1;
    for n = 1:Nstages        
        FA = Min*H0spec.F0; % Sampling frequency of the decimation stage
        HA_fs = HA_fs_v(n);
        HA_Ap = HA_Ap_v(n);
        HA_As = HA_As_v(n);        
        HdA_v{n} = design_filter(FA*1000,H0spec.fp*1000,HA_fs*1000,HA_Ap,HA_As);
        Min = Min / M_v(n);
    end

    fprintf(1,'---------------------------------------\n');
    fprintf(1,'Normalized passband cutoff     : %.4f\n',H0spec.fp/(H0spec.F0/2));
    fprintf(1,'Normalized stopband cutoff     : %.4f\n',H0spec.fs/(H0spec.F0/2));
    fprintf(1,'Single-rate FIR filter length  : %d\n',length(Hd0.Numerator));
    fprintf(1,'Number of stages               : %d\n',Nstages);
    fprintf(1,'Down/upsampling factors        :');
    for n = 1:Nstages, fprintf(1,' %d',M_v(n)); end
    fprintf(1,'\n');
    fprintf(1,'Multirate kernel filter length : %d\n',length(HdK.Numerator));
    fprintf(1,'Antialias/image filter lengths :');
    for n = 1:Nstages, fprintf(1,' %d',length(HdA_v{n}.Numerator)); end
    fprintf(1,'\n');    
    fprintf(1,'---------------------------------------\n');
    
    compare_design_responses(Hd0,M_v,HdK,HdA_v);
    
end

%% Filter design function

function Hd = design_filter(Fs,Fpass,Fstop,Apass,Astop)

    Cf = 0.9; % Correction factor in order to match Dpass with Apass    
    Dpass = Cf * (10^(Apass/40)-1);
    Dstop = 10^(-Astop/20);    
    dens  = 20;               % Density Factor

    % Calculate the filter order
    [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);
    if N < 3, N = 3; end
    
    % Calculate the filter coefficients 
    
    b  = firpm(N, Fo, Ao, W, {dens});
    Hd = dfilt.dffir(b);

end

%% Verification tool

function compare_design_responses(Hd0,M_v,HdK,HdA_v)

    Nstages = length(M_v);
    
    % [1] Impulse input
    Nextra = 2000;
    input = zeros(1,length(Hd0.Numerator)+Nextra);
    input(1) = 1;    
    
    % [2] Simulating single-rate case
    Y0 = filter(Hd0,input);

    % [3] Simulating multirate case    
    YD = input;
    for n = 1:Nstages
        % --- Anti-alias
        YD = filter(HdA_v{n},YD);    
        % --- Downsampling
        YD = YD(1:M_v(n):end);    
    end
    % --- Kernel filtering
    YK = filter(HdK,YD);    
    
    YI = YK;
    for n = Nstages:-1:1
        
        % --- Upsampling
        YI = kron(YI,[1,zeros(1,M_v(n)-1)]);    
        % --- Anti-imaging
        YI = filter(HdA_v{n},YI);            
        YI = M_v(n) * YI;
    end
    
    % [4] Visualization using fvtool    
    N = min(length(Y0),length(YI));
    Y0 = Y0(1:N);
    YI = YI(1:N);        
    R0 = dfilt.dffir;
    R0.Numerator = Y0;
    RI = dfilt.dffir;
    RI.Numerator = YI;    
    fvtool(R0,RI);
    
end

%%