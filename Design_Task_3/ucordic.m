% UCORDIC Unified CORDIC simulator
%
% << Syntax >>
%    [xN,yN,zN,AN] = ucordic(system,mode,x0,y0,z0,N,gaincomp,override,fixptpar,verbose)
%
% << Arguments >>
%    system     'circular', 'hyperbolic', or 'linear'
%    mode       'rotation', or 'vectoring'
%    x0,y0,z0   Initialization parameters : double scalars (see Note-1&2)
%    N          Number of iterations 
%    gaincomp   Flag: Compensate for CORDIC gain? If false then xN and yN
%               have gain in circular and hyperbolic systems.
%    override   Flag: use floating-point values instead of fixed point? 
%    fixptpar   Fixed point word length parameters : structure (see Note-3)
%               Fields:
%               - xy_p, xy_n  Word and fraction lengths used for x and y
%               - z_p, z_n    Same for z 
%    verbose    Flag: show intermediate results?
%    
%    xN,yN,zN   Output values
%    AN         The CORDIC gain
%
% << Description >>
%    The function can be used to simulate various iterative CORDIC 
%    calculations. 
% 
%    Note-1: the implementation does not perform any initial processing to
%            handle large angles. Thus, there are range limits for z0.
%    Note-2: angle values must be given as radiands
%    Note-3: binary point scaling used for all coordinates in this simulation.
%            Computation fails completely if there can be overflow!
%
% << References >>
%    SPS Fall 2014 - Introduction to Part III (2/2) Section 2
%
% << Examples >>
%    % Givens transform
%    fixptpar.xy_p = 32; fixptpar.xy_n = 29;
%    fixptpar.z_p = 32; fixptpar.z_n = 29;
%    x = 1; y = 1.23; phi = 25 * pi / 180;
%    [xd,yd,~,~] = ucordic('circular','rotation',x,y,phi,15,true,false,fixptpar);
%    xt = x * cos(phi) - y * sin(phi); yt = y * cos(phi) + x * sin(phi);
%    fprintf('X: true %.6f cordic %.6f\n',xt,xd);
%    fprintf('Y: true %.6f cordic %.6f\n',yt,yd);

% $Id$

function [xN,yN,zN,AN] = ucordic(system,mode,x0,y0,z0,N,gaincomp,override,fixptpar,verbose)

    if ~exist('verbose','var')
        verbose = false;
    end
    
    % System setup for the given iteration count N
    
    switch system
        case 'circular'
            m = 1;                             % Parameter m = 1
            s = shifts(N,false);               % Amounts of shifts to right
            smult = 2.^(-s);                   % Corresponding multipliers
            e = atan(smult);                   % Rotation angles, arctan
            AN = prod(sqrt(1+2.^(-2 * s)));    % Cordic gain
        case 'hyperbolic'
            m = -1;
            s = shifts(N,true);                % Some steps need repetition
            smult = 2.^(-s);
            e = atanh(smult);
            AN = prod(sqrt(1-2.^(-2 * s)));    
        case 'linear'
            m = 0;
            s = shifts(N,false);
            smult = 2.^(-s);
            e = smult;
            AN = 1;                            % No gain in linear mode
        otherwise
            error('bad system argument (#1)');            
    end
        
    % Simulation    
                
    PREF = fipref;
    if override
        
        % Use double precision to compute ground truth approximation
        % of the intended computation. The result can be used to evaluate
        % the precision of a fixed-point implementation.
        
        PREF.DataTypeOverride = 'TrueDoubles';
        if verbose
            fprintf(1,'Using floating-point override.\n');
        end
    else
        
        % Simulate fixed-point implementation
        
        PREF.DataTypeOverride = 'ForceOff';
        if verbose
            fprintf(1,'Using fixed-point arithmetic.\n');
            fprintf(1,'Format for x and y : s%d.%d\n',fixptpar.xy_p,fixptpar.xy_n);
            fprintf(1,'Format for z       : s%d.%d\n',fixptpar.z_p,fixptpar.z_n);        
        end
    end
        
    % Fixed point modes/attributes to be applied in the computations
    
    FPAxy = fimath('ProductMode','SpecifyPrecision',...
                       'ProductWordLength',fixptpar.xy_p,...
                       'ProductFractionLength',fixptpar.xy_n,...
                       'SumMode','SpecifyPrecision',...
                       'SumWordLength',fixptpar.xy_p,...
                       'SumFractionLength',fixptpar.xy_n,...
                       'RoundMode','fix',...
                       'OverflowMode','wrap');
    FPAz = fimath('ProductMode','SpecifyPrecision',...
                      'ProductWordLength',fixptpar.z_p,...
                      'ProductFractionLength',fixptpar.z_n,...
                      'SumMode','SpecifyPrecision',...
                      'SumWordLength',fixptpar.z_p,...
                      'SumFractionLength',fixptpar.z_n,...
                      'RoundMode','fix',...
                      'OverflowMode','wrap');
                  
    % Conversion to fixed-point representation
    
    x_n = fi(x0,1,fixptpar.xy_p,fixptpar.xy_n);        
    x_n.fimath = FPAxy;
    y_n = fi(y0,1,fixptpar.xy_p,fixptpar.xy_n);        
    y_n.fimath = FPAxy;
    z_n = fi(z0,1,fixptpar.z_p,fixptpar.z_n);        
    z_n.fimath = FPAz;        
    m_fix = fi(m,1,fixptpar.xy_p,fixptpar.xy_n); 
    m_fix.fimath = FPAxy;
    smult_fix = fi(smult,1,fixptpar.xy_p,fixptpar.xy_n); 
    smult_fix.fimath = FPAxy;
    e_fix = fi(e,1,fixptpar.z_p,fixptpar.z_n); 
    e_fix.fimath = FPAz;
        
    % Iterate
    
    if verbose
        fprintf(1,'[Init] X: %.6f Y:%.6f Z:%.6f\n',double(x_n),double(y_n),double(z_n));
    end
    switch mode
          case 'rotation'
              for iter = 1:N
                  x_c = x_n;
                  y_c = y_n;
                  z_c = z_n; 
                  if z_c < 0
                      x_n = x_c + m_fix * y_c * smult_fix(iter);
                      y_n = y_c - x_c * smult_fix(iter);
                      z_n = z_c + e_fix(iter);
                  else
                      x_n = x_c - m_fix * y_c * smult_fix(iter);
                      y_n = y_c + x_c * smult_fix(iter);
                      z_n = z_c - e_fix(iter);
                  end
                  if verbose
                      fprintf(1,'[%d] X: %.6f Y:%.6f Z:%.6f\n',iter,double(x_n),double(y_n),double(z_n));
                  end
              end
              
          case 'vectoring'
              for iter = 1:N
                  x_c = x_n;
                  y_c = y_n;
                  z_c = z_n;
                  if y_c > 0
                      x_n = x_c + m_fix * y_c * smult_fix(iter);
                      y_n = y_c - x_c * smult(iter);
                      z_n = z_c + e(iter);
                  else
                      x_n = x_c - m_fix * y_c * smult_fix(iter);
                      y_n = y_c + x_c * smult_fix(iter);
                      z_n = z_c - e_fix(iter);
                  end
                  if verbose
                      fprintf(1,'[%d] X: %.6f Y:%.6f Z:%.6f\n',iter,double(x_n),double(y_n),double(z_n));
                  end
              end
            otherwise
                error('bad mode argument (#2)');
    end
                
    % Gain compensation (if wanted). Using just fixed-point multiplication under the given format

    zN = double(z_n);
    if gaincomp            
        inv_AN = 1/AN;
        inv_AN_fix = fi(inv_AN,1,fixptpar.xy_p,fixptpar.xy_n);
        xN = double(x_n * inv_AN_fix);
        yN = double(y_n * inv_AN_fix);
        if verbose
            fprintf(1,'[No gain] X: %.6f Y:%.6f Z:%.6f\n',xN,yN,zN);
        end
    else
        xN = double(x_n);
        yN = double(y_n);            
    end        
        
end

%%

function s = shifts(N,repeat)

    if ~repeat
        s = 0:(N-1);
    else
        % Hyperbolic system requires some repetition in order
        % to converge (Andraka, 1998).
        if N <= 4,
            s = 1:N;
        elseif N <= 14
            s = [1:4,4:(N-1)]; % Repeat I = 4
        elseif N <= 42
            s = [1:4,4:13,13:(N-2)]; % Repeat I = 4 and I = 13
        else
            error('too large N'); % Not implemented
        end
    end
    
end

%%
