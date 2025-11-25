% CIC    Tool for simulating Cascaded Integrator Comb (CIC) systems
%
% << Syntax >>
%    sysfun = cic('filter',N,D)
%    sysfun = cic('decim',N,R,M)               % CIC decimator
%    sysfun = cic('interp',N,R,M)              % CIC interpolator
%
% << Arguments >>
%    N        Number of integrator and comb stages (N >= 1)
%    D        Differential delay of comb stage(s) (D >= 1)
%    R        Down/upsampling factor (R >= 2)
%    M        Filter design parameter (typ. either 1 or 2)
%    sysfun   Handle of the system simulation function
%
% << Description >>
%    Tool for simulating various CIC system configurations. 
%    The Fixed-Point Toolbox is needed for using the function.
%
%    The call syntax of the simulation function "sysfun" is
%
%         [output,iplot] = sysfun(input,Nzeropad,p_i,p_o,n)
%
%    where "input" is the vector of input signal samples, "Nzeropad" is
%    the number of zeros appended to the input, and "output" is the vector
%    of system output signal samples, whose length is the length of (input
%    + Nzeropad) x rate change. The second output, "iplot", is a function
%    handle for viewing the signals at intermediate stages: syntax of call
%
%        iplot(phase), where phase = 1 or 2
%
%    The parameters "p_i" and "p_o" denote the input and output word
%    lengths, and "n" is the number of fraction bits. Note that the system
%    does not work correctly, if the setting of these parameters is not
%    correct. 
%
%    For example, the first 1000 samples of CIC filter (N=1, D=10)
%    response to sequence of ones (i.e. step response) can be obtained with the calls
%
%         sysfun = cic('filter',1,10)
%         output = sysfun(ones(1,1000),0,10,14,8)
%
%    where the format of input has been set to s10.8, and
%    output format is the one obtained using Hogenauer's theory, s14.8.
%
%    Other example: CIC decimator (N=2, R=3, M=2)
%
%         sysfun = cic('decim',2,3,2)
%         [output,iplot] = sysfun(ones(1,1000),0,10,16,8)
%         figure(1); iplot(1)
%         figure(2); iplot(2)
%

function sysfun = cic(type,varargin)
    switch type
        case 'filter'
            sysfun = p_cic_filter(varargin{:});
        case 'decim'
            sysfun = p_cic_decimator(varargin{:});
        case 'interp'
            sysfun = p_cic_interpolator(varargin{:});
        otherwise
            error('bad systype argument (#1). Should be ''filter'', ''decim'' or ''interp''');
    end    
end

%%
function sysfun = p_cic_filter(N,D)
    assert(nargin == 2,'two input arguments required: N, D');
    assert(p_isint(N) && p_isint(D),'inputs N and D must be integers');
    assert(N >= 1,'N must be at least 1');
    assert(D >= 1,'D must be at least 1');
    sysfun = @x_sysfun;        
    function [output,iplot] = x_sysfun(input,Nzeropad,p_i,p_o,n,combs_first)
        if nargin < 6, combs_first = false; end
        [input_fi,interm_fi,F] = p_prepare_io(input,Nzeropad,p_i,p_o,n);

        % Due to the linearity of integrators and combs, we
        % simply apply first integrators and then combs. You can
        % also try out alternative order by setting the "combs_first"
        % flag in the call
        %
        % Note: in intro4a.pdf Figure 8(a), the stages are            
        % displayed in mixed order

        valuecollect = p_interm_value_collector(N,combs_first);        
        iplot = valuecollect.plot;
        if ~combs_first
            % Integrator stages first            
            integrator_fifos = p_initialize_fifos(N,1,p_o,n,F);            
            Nticks = length(input_fi);
            valuecollect.count(1,Nticks);
            for tick = 1:Nticks
                integrator_fifos = p_shift_fifos(integrator_fifos);                
                prev_value = input_fi(tick);
                valuecollect.put(1,0,tick,prev_value);
                for stage = 1:N
                    prev_value = prev_value + integrator_fifos{stage}(1);                    
                    integrator_fifos{stage}(end) = prev_value;
                    valuecollect.put(1,stage,tick,prev_value);
                end
                interm_fi(tick) = prev_value;
            end            
            
            % Then combs
            comb_fifos = p_initialize_fifos(N,D,p_o,n,F);            
            output = zeros(1,Nticks);
            valuecollect.count(2,Nticks);
            for tick = 1:Nticks
                comb_fifos = p_shift_fifos(comb_fifos);                
                prev_value = interm_fi(tick);
                valuecollect.put(2,0,tick,prev_value);
                for stage = 1:N
                    comb_fifos{stage}(end) = prev_value;
                    prev_value = prev_value - comb_fifos{stage}(1);                                        
                    valuecollect.put(2,stage,tick,prev_value);
                end
                output(tick) = prev_value;
            end                                    
        else
            % Comb stages first            
            comb_fifos = p_initialize_fifos(N,D,p_o,n,F);            
            Nticks = length(input_fi);
            valuecollect.count(1,Nticks);
            for tick = 1:Nticks
                comb_fifos = p_shift_fifos(comb_fifos);                
                prev_value = input_fi(tick);
                valuecollect.put(1,0,tick,prev_value);
                for stage = 1:N
                    comb_fifos{stage}(end) = prev_value;
                    prev_value = prev_value - comb_fifos{stage}(1);                                        
                    valuecollect.put(1,stage,tick,prev_value);
                end
                interm_fi(tick) = prev_value;
            end                        

            % Then integrators
            integrator_fifos = p_initialize_fifos(N,1,p_o,n,F);            
            output = zeros(1,Nticks);
            valuecollect.count(2,Nticks);
            for tick = 1:Nticks
                integrator_fifos = p_shift_fifos(integrator_fifos);                
                prev_value = interm_fi(tick);
                valuecollect.put(2,0,tick,prev_value);
                for stage = 1:N
                    prev_value = prev_value + integrator_fifos{stage}(1);                    
                    integrator_fifos{stage}(end) = prev_value;
                    valuecollect.put(2,stage,tick,prev_value);
                end
                output(tick) = prev_value;
            end            
        end
    end
end

%%
function sysfun = p_cic_decimator(N,R,M)
    assert(nargin == 3,'three input arguments required: N, R, M');
    assert(p_isint(N) && p_isint(R) && p_isint(M),'inputs N, R, and M must be integers');
    assert(N >= 1,'N must be at least 1');
    assert(R >= 2,'R must be at least 2');
    assert(M >= 1,'M must be at least 1');
    sysfun = @x_sysfun;
    function [output,iplot] = x_sysfun(input,Nzeropad,p_i,p_o,n)
        [input_fi,interm_fi,F] = p_prepare_io(input,Nzeropad,p_i,p_o,n);
        valuecollect = p_interm_value_collector(N,false);        
        iplot = valuecollect.plot;

        % Integrator stages first            
        integrator_fifos = p_initialize_fifos(N,1,p_o,n,F);            
        Nticks = length(input_fi);
        valuecollect.count(1,Nticks);
        for tick = 1:Nticks
            integrator_fifos = p_shift_fifos(integrator_fifos);                
            prev_value = input_fi(tick);
            valuecollect.put(1,0,tick,prev_value);
            for stage = 1:N
                prev_value = prev_value + integrator_fifos{stage}(1);                    
                integrator_fifos{stage}(end) = prev_value;
                valuecollect.put(1,stage,tick,prev_value);
            end
            interm_fi(tick) = prev_value;
        end            

        % Downsampling
        interm_fi = downsample(interm_fi,R);

        % Then combs
        comb_fifos = p_initialize_fifos(N,M,p_o,n,F);            
        Nticks = length(interm_fi);
        output = zeros(1,Nticks);
        valuecollect.count(2,Nticks);
        for tick = 1:Nticks
            comb_fifos = p_shift_fifos(comb_fifos);                
            prev_value = interm_fi(tick);
            valuecollect.put(2,0,tick,prev_value);
            for stage = 1:N
                comb_fifos{stage}(end) = prev_value;
                prev_value = prev_value - comb_fifos{stage}(1);                                        
                valuecollect.put(2,stage,tick,prev_value);
            end
            output(tick) = prev_value;
        end                        
    end
end

%%
function sysfun = p_cic_interpolator(N,R,M)
    assert(nargin == 3,'three input arguments required: N, R, M');
    assert(p_isint(N) && p_isint(R) && p_isint(M),'inputs N, R, and M must be integers');
    assert(N >= 1,'N must be at least 1');
    assert(R >= 2,'R must be at least 2');
    assert(M >= 1,'M must be at least 1');
    sysfun = @x_sysfun;
    function [output,iplot] = x_sysfun(input,Nzeropad,p_i,p_o,n)
        [input_fi,interm_fi,F] = p_prepare_io(input,Nzeropad,p_i,p_o,n);
        valuecollect = p_interm_value_collector(N,true);        
        iplot = valuecollect.plot;

        % Comb stages first            
        comb_fifos = p_initialize_fifos(N,M,p_o,n,F);            
        Nticks = length(input_fi);
        valuecollect.count(1,Nticks);
        for tick = 1:Nticks
            comb_fifos = p_shift_fifos(comb_fifos);                
            prev_value = input_fi(tick);
            valuecollect.put(1,0,tick,prev_value);
            for stage = 1:N
                comb_fifos{stage}(end) = prev_value;
                prev_value = prev_value - comb_fifos{stage}(1);                                        
                valuecollect.put(1,stage,tick,prev_value);
            end
            interm_fi(tick) = prev_value;
        end                        

        % Upsampling            
        interm_fi = upsample(interm_fi,R);

        % Then integrators
        integrator_fifos = p_initialize_fifos(N,1,p_o,n,F);            
        Nticks = length(interm_fi);
        output = zeros(1,Nticks);
        valuecollect.count(2,Nticks);

        for tick = 1:Nticks
            integrator_fifos = p_shift_fifos(integrator_fifos);                
            prev_value = interm_fi(tick);
            valuecollect.put(2,0,tick,prev_value);
            for stage = 1:N
                prev_value = prev_value + integrator_fifos{stage}(1);                    
                integrator_fifos{stage}(end) = prev_value;
                valuecollect.put(2,stage,tick,prev_value);
            end
            output(tick) = prev_value;
        end                        
    end
end

%%

function [input_fi,interm_fi,F] = p_prepare_io(input,Nzeropad,p_i,p_o,n)
    F = fimath('SumMode','SpecifyPrecision','SumWordLength',p_o,'SumFractionLength',n,'OverflowAction','Wrap');
    input = [reshape(input,[1,numel(input)]),zeros(1,Nzeropad)];
    input_fi = fi(input,1,p_i,n,F); 
    interm_fi = fi(zeros(size(input)),1,p_o,n,F);
end

%%

function ok = p_isint(val)
    ok = mod(val,1) == 0;
end

%%

function [fifos] = p_initialize_fifos(N,D,p_o,n,F)
    fifos = cell(1,N);
    for ind = 1:N
        fifos{ind} = fi(zeros(1,D+1),1,p_o,n,F);
    end        
end

function fifos = p_shift_fifos(fifos)
    N = length(fifos);    
    for ind = 1:N
        fifos{ind}(1:(end-1)) = fifos{ind}(2:end);
        fifos{ind}(end) = 0;
    end
end

%%

function h = p_interm_value_collector(Nstages,combs_first)

    phase1 = [];
    phase2 = [];
    
    h.count = @x_count;
    h.put = @x_put;
    h.plot = @x_plot;
    
    function x_count(phase,Nticks)
        if phase == 1
            phase1 = nan(Nstages+1,Nticks);
        else
            phase2 = nan(Nstages+1,Nticks);
        end        
    end

    function x_put(phase,stage,tick,value)
        if phase == 1
            phase1(stage+1,tick) = value;
        else
            phase2(stage+1,tick) = value;
        end
    end

    function x_plot(phase)
        titles = cell(1,Nstages+1);
        if phase == 1
            if combs_first
                titles{1} = 'Comb stages: input';
            else
                titles{1} = 'Integrator stages: input';
            end
            ph = phase1;
        elseif phase == 2
            if combs_first
                titles{1} = 'Integrator stages: input';
            else
                titles{1} = 'Comb stages: input';
            end        
            ph = phase2;
        end
        for stage = 1:Nstages
            titles{stage+1} = sprintf('Stage %d output',stage);
        end
        if phase == 1
            titles{1} = sprintf('%s (= system input)',titles{1});
        end        
        if phase == 2
            titles{Nstages+1} = sprintf('%s (= system output)',titles{Nstages+1});
        end
        for stage = 0:Nstages
            subplot(Nstages+1,1,stage+1);
            plot(ph(stage+1,:));
            set(gca,'xlim',[0,size(ph,2)+1],'ylim',[min(ph(stage+1,:))-0.1,max(ph(stage+1,:))+0.1]);
            title(titles{stage+1});
        end
    end

end
