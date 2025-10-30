% NN_FCRELUFC_SIMext      Neural network: FC+RELU+FC simulator (extended version)
%
% << Syntax >>
%    h = NN_FCRELUFC_SIM(flparams,hwprops)
%
% << Arguments >>
%    flparams      Floating-point NN parameters, struct with fields
%                  - weights_0 : Layer 0 weights
%                  - bias_0 : Layer 0 bias
%                  - weights_1 : Layer 1 weights
%                  - bias_1 : Layer 1 bias
%    hwprops       Properties of the NN implementation hardware, struct with fields
%                  - NbitsIO   Input/output word length. Same word length available for
%                              storing weights. For biases, doubled length is used.
%                  - NbitsW    Number of bits available for weights
%                  - Nguard    Number of guard bits in the accumulator, that is,
%                              MAC accumulator word length = NbitsIO + NbitsW + Nguard
%                  - postQsas  Shift-add count for post-quantization of MAC output.
%                              
%                    The idea is that before extraction of bits, the output is 
%                    multiplied using shift-add logic. Also, we assume that 
%                    most-significant bits are extracted by HW after multiplication.
%                    The goal of multiplication is to increase utilization of 
%                    available range provided by "NbitsIO" bits. 
%
% << Description >>
%    Simulator for a NN, which consists of input layer, fully connected
%    hidden layer, RelU activation layer, and fully connected output layer.
%
%    HW includes MAC units, which use asymmetrically quantized weights as input. 
%    MAC units contain also CORDIC based multipliers for adjusting accumulator outputs 
%    for I/O word length. 
%
%    Function
%
%        [calib_flo,calib_int] = h.calibrate(calib_input,max_in,Nbits_flo)
%
%    adjusts the HW for operating with the specified input. "calib_flo"
%    provides input and output of the original floating-point filter for
%    the calib_input quantized to I/O word length. "calib_int" provides
%    corresponding input and output for the quantized NN. If "Nbits_flo"
%    is given, then "calib_flo.input" is quantized using this word length,
%    instead of using I/O word length.
%
%    Function
%
%        [data_flo,data_int] = h.run(data_input)
%
%    executes calibrated NN for given data. "data_flo" provides input and 
%    output of the original floating-point filter for the data_input quantized 
%    to I/O word length.

function h = NN_FCRELUFC_SIM(flparams,hwprops)

    [input_len,flparams] = check_flparams(flparams); % input_len is the length of the NN input vector
    hwprops.NbitsACCU = hwprops.NbitsIO + hwprops.NbitsW + hwprops.Nguard;

    max_input = [];
    nnparams = [];
    NNoutput_scale = [];

    h.calibrate = @x_calibrate;
    h.run = @x_run;
    h.evaluate = @x_evaluate; % Evaluate performance for input data

    h.nnparams = @x_nnparams;
    h.visualize_ranges = @x_visualize_ranges;

    function x_visualize_ranges
        subplot(211);
        plot(nnparams.weights_0(:),flparams.weights_0(:),'.');
        set(gca,'xlim',[-0.3,2^hwprops.NbitsW-1+0.3]);
        set(gca,'xgrid','on','ygrid','on');
        title('Layer 0 weight encoding')

        subplot(212);
        plot(nnparams.weights_1(:),flparams.weights_1(:),'.');
        set(gca,'xlim',[-0.3,2^hwprops.NbitsW-1+0.3]);
        set(gca,'xgrid','on','ygrid','on');
        title('Layer 1 weight encoding')
        
    end

    % Calibrate
    %
    %     calib_input  Good set of representative samples, values must be
    %                  within the range [-max_in, +max_in], where max_in
    %                  corresponds to the integer 2^(NbitsIO-1)-1. 
    %
    %     Output "calib_flo.input" provides exact floating point values,
    %     which correspond to integers obtained from "calib_input" conversion
    %     Output "calib_flo.output" provides NN floating point output for
    %     that input.

    function [calib_flo,calib_int] = x_calibrate(calib_input,max_in,Nbits_flo)

        fprintf('* HW properties *\n')
        fprintf('I/O word length: %d (used for all NN layers)\n',hwprops.NbitsIO);
        fprintf('Weight encoding: %d bits, asymmetric\n',hwprops.NbitsW);
        fprintf('ACCU word length: %d (%d guard bits)\n',hwprops.NbitsACCU,hwprops.Nguard);
        fprintf('ACCU output mapping, gain shift-add count: %d\n',hwprops.postQsas);
        fprintf('\n');

        if nargin < 3
            fprintf('No input quantization for floating-point NN.\n');
        else
            if Nbits_flo < hwprops.NbitsIO
                error('Nbits_flo (3rd argument) must be at least hwprops.NbitsIO\n');
            end
            fprintf('Input quantization for floating-point NN uses %d bits.\n',Nbits_flo);
        end
        
        fprintf('* Calibration *\n')
        max_input = max_in;
        assert(size(calib_input,1) == input_len);
        assert(all(abs(calib_input(:)) <= max_input));

        % Map calibration data to integer input
        maxint = 2^(hwprops.NbitsIO-1)-1;
        nnparams.input.s = max_input/maxint;
        calib_int.input = round((1/nnparams.input.s)*calib_input);

        if nargin < 3
            calib_flo.input = single(nnparams.input.s * calib_int.input);
        else
            maxint_flo = 2^(Nbits_flo-1)-1;
            flo_s = max_input/maxint_flo;
            calib_flo.input = single(flo_s * round((1/flo_s)*calib_input));
        end

        fprintf('DIFF provide information about the quality of results at various stages of NN computation.\n\n')
        
        i_compare('input',calib_int.input,calib_flo.input,nnparams.input.s); % Perfect match
        i_note('no difference because floating-point input is tuned to match with integer input (for fair comparison)')

        Ncalib_data = size(calib_input,2);

        l0_out = flparams.weights_0 * calib_flo.input + repmat(flparams.bias_0,1,Ncalib_data);
        relu_out = max(l0_out,0);
        calib_flo.output = flparams.weights_1 * relu_out + repmat(flparams.bias_1,1,Ncalib_data);

        calib_flo.range_layer0 = [min(l0_out(:)),max(l0_out(:))];
        calib_flo.range_relu = [0,max(relu_out(:))];
        calib_flo.range_output = [min(calib_flo.output(:)),max(calib_flo.output(:))];
        
        % Compute asymmetric weights for fc layer 0
        subplot(211);
        [s,zeropt,w_uint] = search_asymmetric_parameters(hwprops.NbitsW,flparams.weights_0,'minmse');
        title('FC layer 0 scale search');
        nnparams.encode_0.s = s;
        nnparams.encode_0.zeropt = zeropt;
        nnparams.encode_0.Nbits = hwprops.NbitsW;        
        nnparams.weights_0 = w_uint;
        fprintf('Asymmetric coding for layer 0 weights: scale %.6f zeropt %d Nbits %d - uint/range %d..%d\n',s,zeropt,hwprops.NbitsW,min(w_uint(:)),max(w_uint(:)));

        % Compute integers for bias (symmetric encoding used)
        accu0_scale = nnparams.input.s * nnparams.encode_0.s; 
        nnparams.bias_0 = round((1/accu0_scale) * flparams.bias_0);
        fprintf('Symmetric coding for layer 0 bias: scale %.6f - int/range %d..%d\n',accu0_scale,min(nnparams.bias_0),max(nnparams.bias_0));

        % Simulate layer_0 MAC to get accumulator values
        % This data is needed to adjust quantization of layer_0 output
        [accuvals,Nsaturate] = simulate_fclayer(calib_int.input,nnparams.weights_0,nnparams.encode_0.zeropt,nnparams.bias_0,hwprops.NbitsACCU);
        if Nsaturate > 0, warning('layer 0: accumulator saturates for %d (%.2f%%) inputs',Nsaturate,100*Nsaturate/numel(accuvals)); end

        i_compare('layer0-before-Q (i.e. accu)',accuvals,l0_out,accu0_scale);

        % Define post-quantization setup for layer 0
        nnparams.postQ_0 = define_post_quantization(accuvals,hwprops.NbitsACCU,hwprops.NbitsIO,hwprops.postQsas);

        % Produce layer_0 output
        [postQvals,scalingQ0] = apply_postQ(accuvals,nnparams.postQ_0);
        nnparams.postQ_0.s = 1 / scalingQ0;
        layer0output_scale = accu0_scale * nnparams.postQ_0.s;

        i_compare('layer0-after-Q (i.e. output)',postQvals,l0_out,layer0output_scale);
        i_note('Compared to previous difference, increase in difference (less precision)');
        
        % Produce RelU output
        RelUout = max(0,postQvals);
        RelUoutput_scale = layer0output_scale;

        i_compare('relu-out',RelUout,relu_out,RelUoutput_scale);
        
        % Compute asymmetric weights for fc layer 1
        subplot(212);
        [s,zeropt,w_uint] = search_asymmetric_parameters(hwprops.NbitsW,flparams.weights_1,'minmse');
        title('FC layer 1 scale search');
        nnparams.encode_1.s = s;
        nnparams.encode_1.zeropt = zeropt;
        nnparams.encode_1.Nbits = hwprops.NbitsW;        
        nnparams.weights_1 = w_uint;
        fprintf('Asymmetric coding for layer 1 weights: scale %.6f zeropt %d Nbits %d - uint/range %d..%d\n',s,zeropt,hwprops.NbitsW,min(w_uint(:)),max(w_uint(:)));
        
        % Compute integers for bias (symmetric encoding used)
        accu1_scale = RelUoutput_scale * nnparams.encode_1.s;
        nnparams.bias_1 = round((1/accu1_scale) * flparams.bias_1);
        fprintf('Symmetric coding for layer 1 bias: scale %.6f - int/range %d..%d\n',accu1_scale,min(nnparams.bias_1),max(nnparams.bias_1));

        % Simulate layer_1 MAC to get accumulator values
        [accuvals,Nsaturate] = simulate_fclayer(RelUout,nnparams.weights_1,nnparams.encode_1.zeropt,nnparams.bias_1,hwprops.NbitsACCU);
        if Nsaturate > 1, warning('layer 1: accumulator saturates for %d (%.2f%%) inputs',Nsaturate,100*Nsaturate/numel(accuvals)); end

        i_compare('layer1-before-Q (i.e. accu)',accuvals,calib_flo.output,accu1_scale);
        
        % Define post-quantization setup for layer 0
        nnparams.postQ_1 = define_post_quantization(accuvals,hwprops.NbitsACCU,hwprops.NbitsIO,hwprops.postQsas);

        % Produce layer_0 output
        [calib_int.output,scalingQ1] = apply_postQ(accuvals,nnparams.postQ_1);
        nnparams.postQ_1.s = 1 / scalingQ1;
        NNoutput_scale = accu1_scale * nnparams.postQ_1.s;

        i_compare('layer1-after-Q (i.e. NN output)',calib_int.output,calib_flo.output,NNoutput_scale);
        i_note('Compared to previous difference, increase in difference (less precision)');
        
    end

    function [data_flo,data_int] = x_run(data_input)
        if isempty(nnparams), error('not calibrated yet'); end

        assert(size(data_input,1) == input_len);

        % Clipping input if going out of range
        I = data_input < -max_input;
        data_input(I) = -max_input;
        I = data_input > max_input;
        data_input(I) = max_input;

        % Map data to integer input
        data_int.input = round((1/nnparams.input.s)*data_input);
        data_flo.input = single(nnparams.input.s * data_int.input);

        % i_compare('input',data_int.input,data_flo.input,nnparams.input.s); % Perfect match
        
        Ndata = size(data_input,2);

        l0_out = flparams.weights_0 * data_flo.input + repmat(flparams.bias_0,1,Ndata);
        relu_out = max(l0_out,0);
        data_flo.output = flparams.weights_1 * relu_out + repmat(flparams.bias_1,1,Ndata);

        [accuvals,Nsaturate] = simulate_fclayer(data_int.input,nnparams.weights_0,nnparams.encode_0.zeropt,nnparams.bias_0,hwprops.NbitsACCU);
        if Nsaturate > 0, warning('layer 0: accumulator saturates for %d (%.2f%%) inputs',Nsaturate,100*Nsaturate/numel(accuvals)); end

        [postQvals,~] = apply_postQ(accuvals,nnparams.postQ_0);

        RelUout = max(0,postQvals);

        [accuvals,Nsaturate] = simulate_fclayer(RelUout,nnparams.weights_1,nnparams.encode_1.zeropt,nnparams.bias_1,hwprops.NbitsACCU);
        if Nsaturate > 1, warning('layer 1: accumulator saturates for %d (%.2f%%) inputs',Nsaturate,100*Nsaturate/numel(accuvals)); end

        [data_int.output,~] = apply_postQ(accuvals,nnparams.postQ_1);

        i_compare('layer1-after-Q (i.e. NN output)',data_int.output,data_flo.output,NNoutput_scale);
        
    end

    function v = x_evaluate(data_input)
        [data_flo,data_int] = x_run(data_input);
        v.diff = NNoutput_scale * data_int.output(:) - data_flo.output(:);
        v.output_gt = data_flo.output(:);
        v.output_int = NNoutput_scale * data_int.output(:);
        v.rmse = sqrt(mean(v.diff(:).^2));
        v.mae = mean(abs(v.diff(:)));
    end

    function v = x_nnparams
        if isempty(nnparams), error('not calibrated yet'); end
        v = nnparams;
    end

    function [output,Nsaturate] = simulate_fclayer(input,asym_weights,zeropt,bias,NbitsACCU)
        % Prepare asymmetric weights
        weights = asym_weights - zeropt;
        output = weights * input + repmat(bias,1,size(input,2));
        [output,Nsaturate] = clamp_int(output,NbitsACCU);
    end

    function [x,Nsaturate] = clamp_int(x,Nbits)
        PMAX = 2^(Nbits-1)-1;
        NMAX = -2^(Nbits-1);
        IclipP = (x > PMAX);
        IclipN = (x < NMAX);
        x(IclipP) = PMAX;
        x(IclipN) = NMAX;
        Nsaturate = sum(IclipP(:)) + sum(IclipN(:));       
    end

    function [output_values,scaling] = apply_postQ(accuvals,postQ)
        scaling = postQ.mplier * 2.^(postQ.NbitsOut - postQ.NbitsIn);
        output_values = round(accuvals * scaling); % floor and 2.^x simulate extraction of most significant bits
    end

    function i_compare(txt,int_vals,flo_vals,s)
        diff = s * int_vals(:) - flo_vals(:);
        fprintf('# DIFF at ''%s'': max(abs) %.6f avg %.6f stdev %.6f flo/range %.6f...%.6f\n',txt,max(abs(diff)),mean(diff),std(diff),min(flo_vals(:)),max(flo_vals(:)));
    end

    function i_note(txt)
        fprintf('# NOTE: %s\n',txt);
    end
end

function [input_len,flparams] = check_flparams(flparams)

    assert(all(isfield(flparams,{'weights_0','bias_0','weights_1','bias_1'})));
    [N1,M1] = size(flparams.weights_0);
    [N2,M2] = size(flparams.weights_1);
    assert(N1 == M2);
    assert(all(size(flparams.bias_0) == [N1,1]));
    assert(all(size(flparams.bias_1) == [N2,1]));
    input_len = M1;

    % Conversion to single precision used in floating point simulation
    flparams.weights_0 = single(flparams.weights_0);
    flparams.bias_0 = single(flparams.bias_0);
    flparams.weights_1 = single(flparams.weights_1);
    flparams.bias_1 = single(flparams.bias_1);
    
end

function postQ = define_post_quantization(output,NbitsIn,NbitsOut,postQsas)

    postQ.NbitsIn = NbitsIn;
    postQ.NbitsOut = NbitsOut;
    if postQsas == 0
        postQ.shifts = [];
        postQ.mplier = 1;
        return;
    end

    max_out = max(abs(output(:)));
    target_out = 2^(NbitsIn-1)-1;
    mplier = target_out / max_out;

    % Find best sequence of shifts for approximating "mplier"
    % using the available shift/add resources
    fprintf('Multiplier to be approximated: %.6f\n',mplier);

    postQ.shifts = zeros(1,postQsas);
    mplierx = mplier;
    for n = 1:postQsas
        postQ.shifts(n) = floor(log2(mplierx));
        mplierx = mplierx-2^postQ.shifts(n);
    end
    postQ.mplier = sum(2.^postQ.shifts);
    fprintf('Shifts:%s',sprintf(' %d',postQ.shifts));
    fprintf(' implement %.6f\n',postQ.mplier);

end
