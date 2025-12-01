% CYCLIC_CONV Computes cyclic convolution for a signal section
%
% << Syntax >>
% y_section = cyclic_conv(u_section,h,N,domain)
%
% << Arguments >>
% u_section Section of the input signal (size 1 x Nu)
% h Filter h (size 1 x Nh)
% N The cyclic convolution length (number of points in DFT)
% - this is the length of the output "y"
% domain 'time' (default) or 'frequency' (should give the same output)
% y Convolution result (size 1 x N)
%
% << Description >>
% Computation of cyclic convolution. Can be used to compute
% intermediate values required by sectioning (overlap-add or overlap-save).
% The procedure is:
%
% [1] extract subsequences from the input "u" according to the sectioning rules
% [2] compute cyclic convolution of each subsequence using this function
% [3] combine convolution results according to the sectioning rules
%
% See Introduction Sec. 2.1.3 + Figures 4 and 5.
%
% Note: if sectioning works then the output should correspond to
% "c = conv(h,u)" in the valid part!
 
% $Id: cyclic_conv.m,v 1.2 2013/11/08 07:12:43 psangi Exp $
 
function y = cyclic_conv(u_section,h,N,domain)
 
  if ~exist('domain','var'), domain = 'time'; end
 
  Nu = size(u_section,2);
  Nh = size(h,2);
 
  if (N < Nu) || (N < Nh)
  error('N must greater than or equal to length of section/filter');
  end
 
  % Zero padding
 
  u_section_p = [u_section, zeros(1,N-Nu)];
  h_p = [h, zeros(1,N-Nh)];
 
  % Cyclic convolution: in practice, frequency domain is used, but
  % also time domain computation is
 
  switch domain
 
  case 'frequency'
  % Computation in frequency domain using Matlab's "fft" and "ifft"
  % Corresponds to Equation 5 in Introduction.
  H = fft(h_p);
  U = fft(u_section_p);
  Y = H .* U;
  y = ifft(Y);
 
  case 'time'
  % Computation in time-domain using Matlab's "conv"
  u_section_x = u_section_p([(N-Nh+2):N,1:N]);
  y = conv(u_section_x,h,'valid');
 
  otherwise
  error('bad domain argument');
  end
 
end 
