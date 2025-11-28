H0spec.F0 = 110;
H0spec.fp = 2.8;
H0spec.fs = 3.9;
H0spec.Ap = 0.1;
H0spec.As = 60; 
M_v = [2, 4];
HK_Ap = 0.1;          
HK_As = 60;           
HA_fs_v = [51.1, 9.85];  
HA_Ap_v = [0.002, 0.002];
HA_As_v = [61, 61];

multirateN(H0spec, M_v, HK_Ap, HK_As, HA_fs_v, HA_Ap_v, HA_As_v)

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