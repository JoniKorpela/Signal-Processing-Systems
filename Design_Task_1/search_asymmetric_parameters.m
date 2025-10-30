function [s,zeropt,w_uint] = search_asymmetric_parameters(Nbits,weights,criterion)

    switch criterion
        case 'minmaxlim'
            [s,zeropt] = search_minmaxlim(Nbits,min(weights(:)),max(weights(:)));
        case 'minmse'
            [s,zeropt] = search_minmse(Nbits,weights);            
    end

    [w_uint] = compute_qweights(weights,Nbits,s,zeropt);
    
end

function [s,zeropt] = search_minmaxlim(Nbits,wmin,wmax)

    % Target:
    %     -sz = wmin
    %     s(maxint-z) = wmax
    % where 
    %     maxint = 2^Nbits-1

    maxint = 2^Nbits-1;

    s = (wmax - wmin) / maxint;

    zeropt = round(-wmin/s); % Must be integer

    s = -wmin/zeropt;

end

function [s,zeropt] = search_minmse(Nbits,weights)

    % Target:
    %     -sz = wmin
    %     s(maxint-z) = wmax
    % where 
    %     maxint = 2^Nbits-1

    maxint = 2^Nbits-1;

    wmin = min(weights(:));
    wmax = max(weights(:));

    s0 = (wmax - wmin) / maxint;
    zeropt = round(-wmin/s0); 

    % Use this as the zero point and
    % search for scaling, that minimizes MSE
    % Scalings are selected from a grid controlled
    % by the user. 

    s0 = -wmin/zeropt; % One option to consider (sa
    
    s_mpliers = 0.9:0.01:2.0;
    ss = s0 * s_mpliers; % Tested scale values

    N = numel(s_mpliers);
    mse = zeros(1,N);    
    for ind = 1:N
        s = ss(ind);
        [~,w_approx] = compute_qweights(weights,Nbits,s,zeropt);
        mse(ind) = mean((weights(:) - w_approx(:)).^2);
    end
    [minmse,I] = min(mse);
    s = ss(I(1));
    plot(ss,mse,'.-');
    hold on 
    plot(s,minmse,'rs')
    hold off
    xlabel('scale'); ylabel('mse');
end

%%%%%%%%%%%%%%%%%%
%
% Parameters
% weights  The original floating point weights
% Nbits    Number of bits used in quantization
% s        Scaling factor
% zeropt   Zero point
%
% Returns
% w_uint   Integer quantized weights. Note that you are supposed only to
%          simulate unsigned integer arithmetics, the data type of this 
%          variable must be double.
% w_approx Dequantized coefficients.
%
%%%%%%%%%%%%%%%%%%
function [w_uint,w_approx] = compute_qweights(weights,Nbits,s,zeropt)
    [w_uint, w_approx] = quantizew(weights, Nbits, s, zeropt);
end
