% << Syntax >>
%    h = FIR_MAC_Simulator(qxp,qap,qpp,gbits,qop)
%
% << Arguments >>
%    qxp    Quantization parameters of the filter input (ADC) : [p n f], f optional (default 1)
%    qap    Quantization parameters of the filter coefficients : [p n f], f optional (default 1)
%    qpp    Product word and fraction length : [p n] or [] = use full precision
%    gbits  Accumulator, the number of quard bits
%    qop    The filter output format (post-quantization): [p n]        
%    h      Structure of access function handles
%
% << Description >>
%    The access functions are
%
%        sigall = h.simulate(mode,x,a)
%
%    where "mode" is an index, which selects the simulation mode, "a"
%    contains the filter coefficients, and "x" is the sampled signal.
%    
%    --------------------------------------------------------------------------
%    Value of "mode"               0 1 2 3 4 5 6 7 8 9   Controls
%    --------------------------------------------------------------------------
%    Fixed-point MAC               1 0 0 0 0 0 0 0 0 0   "qpp" and "gbits"
%    Coefficient quantization      1 0 1 1 1 1 0 0 0 0   "qap" 
%    ADC conversion                1 0 0 1 1 0 0 1 1 0   "qxp" 
%    Post-quantization (rounding)  1 0 0 0 1 1 0 0 1 1   "qop"
%    --------------------------------------------------------------------------
%
%    The idea is to allow studying the effect of various quantization
%    operations to the performance of the filter. Floating-point MAC
%    computations are based on double override option. 
%
%       h.

% $Id$

function h = FIR_MAC_Simulator(qxp,qap,qpp,gbits,qop)
    
    if isempty(qpp)
        qpp = qxp(1:2) + qap(1:2); % Full precision : based on the "law of conservation of bits"
    end
    
    % =======================================================================================
    % Fixed point modes/attributes to be applied in the computations
    FPA = fimath('ProductMode','SpecifyPrecision',...
                 'ProductWordLength',qpp(1),...
                 'ProductFractionLength',qpp(2),...
                 'SumMode','SpecifyPrecision',...
                 'SumWordLength',qpp(1)+gbits,...
                 'SumFractionLength',qpp(2),...
                 'RoundMode','nearest',...
                 'OverflowMode','saturate');
             
    signals = [];
        
    h.simulate = @x_simulate;
    h.computeDiff = @x_computeDiff;
    
    %% Computes differences of the current simulation and reference simulation
    

    function diff = x_computeDiff(ref_sigall)
        diff.adc_out = signals.adc_out - ref_sigall.adc_out;
        diff.mac_out = signals.mac_out - ref_sigall.mac_out;
        diff.out = signals.out - ref_sigall.out;        
    end

    %% Performs simulation in the mode selected
    
    function sigall = x_simulate(mode,x,a)
        
        [modeName,fixptMAC,coeffQ,ADC,postQ] = parseMode(mode);
        N = length(x); % Signal length
        L = length(a); % Number of coefficients
        
        fprintf(1,'FIR/MAC Simulation\n');
        fprintf(1,'Execution mode: %s\n',modeName);
        fprintf(1,'Signal length : %d\n',N);
        fprintf(1,'Filter taps   : %d\n',L);        

        PREF = fipref;
        if ~fixptMAC
            PREF.DataTypeOverride = 'TrueDoubles';
        else
            PREF.DataTypeOverride = 'ForceOff';
        end    
                        
        %% =======================================================================================
        % Coefficient quantization
        
        if coeffQ
            if length(qap) < 3, qap(3) = 1; end
            qa = fi(a / qap(3),1,qap(1),qap(2));
            qa.fimath = FPA;
        else
            qa = a;
        end 
        qa = qa(:);
        
        signals.coeff = a(:);
        signals.coeff_quant = double(qa);
        signals.coeff_error = signals.coeff_quant - signals.coeff;
        
        
        signals.in = x(:);
        
        %% =======================================================================================
        % Analog-to-Digital conversion
        % Note that in the floating-point simulation, this must be done without override.
        % For this reason, "quantizer" used here. Maybe just "fi" could be used

        if ADC
            if length(qxp) < 3, qxp(3) = 1; end
            qx = quantizer('Mode','fixed','Roundmode','nearest',...
                           'Overflowmode','saturate','Format',qxp(1:2));
            x = x / qxp(3);
            x = quantize(qx,x);
            qx = fi(x,1,qxp(1),qxp(2));
            qx.fimath = FPA;
        else
            qx = x; % No ADC quantization
        end
        qx = qx(:);
        
        signals.adc_out = double(qx);
        signals.adc_error = signals.adc_out - signals.in;
        
        if fixptMAC
            mac_out = fi(zeros(1,N),1,qpp(1)+gbits,qpp(2));
            y.fimath = FPA;
        else
            mac_out = zeros(1,N);
        end
        
        %% =======================================================================================
        % Reserve space for output and delay line
        
        if postQ
            y = fi(zeros(1,N),1,qop(1),qop(2));
            y.fimath = FPA;
        else
            y = zeros(1,N);
        end

        if ADC
            delay_line = fi(zeros(1,L),1,qxp(1),qxp(2));
            delay_line.fimath = FPA;
        else
            delay_line = zeros(1,L);
        end
        
        %% =======================================================================================
        % Processing of the data sequence

        for k = 1:N
            delay_line = [qx(k), delay_line(1:(L-1))];
            xout = delay_line * qa;
            mac_out(k) = xout; % Accumulator output
            y(k) = mac_out(k); % Rounding occurs here
        end        
        
        PREF.DataTypeOverride = 'ForceOff';    
        
        signals.mac_out = double(mac_out);
        signals.out = double(y);
        signals.trunc_error = signals.out - signals.mac_out;

        sigall = signals;
        
    end

    
end

%%

function [modeName,fixptMAC,coeffQ,ADC,postQ] = parseMode(mode)

    fixptMAC_v = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    coeffQ_v   = [1, 0, 1, 1, 1, 1, 0, 0, 0, 0];
    ADC_v      = [1, 0, 0, 1, 1, 0, 0, 1, 1, 0];
    postQ_v    = [1, 0, 0, 0, 1, 1, 0, 0, 1, 1];
    
    index = mode+1;    
    fixptMAC = fixptMAC_v(index);
    coeffQ = coeffQ_v(index);
    ADC = ADC_v(index);
    postQ = postQ_v(index);

    if fixptMAC
        modeName = 'fixed-point';
    else
        modeName = 'override';
        if coeffQ, modeName = [modeName '+coeffQ']; end
        if ADC, modeName = [modeName '+ADC']; end
        if postQ, modeName = [modeName '+postQ']; end
    end
    
end

%% EOF
