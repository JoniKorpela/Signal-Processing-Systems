% EXSIM_RUN Simulate the exercise model
%
% << Syntax >>
%    exsim_run(L,D,step_size)
%
% << Arguments >>
%    L          LMS filter length
%    D          The delay of the "Integer Delay" block
%    step_size  The step size of the LMS filter (mu')
%
% << Description >>
%    An example of how to start simulation from Matlab. 
%
%    It is convenient to use this kind of function to run simulations
%    if you must try out several possible parameter combinations. 
%
%    NOTE-1: Modify the model first so that it contains references
%            to the variables L, D, and step_size at right places (see
%            design task handout).
%
%    NOTE-2: Function allows systematic testing - call it in a loop which
%            varies the input parameter values.
%
% << Examples >>
%    exsim_run(300,150,0.001)

function [n_errors, pass, cvrg, MSE] = exsim_run(L,D,step_size) %#ok<INUSD>
        
    stop_time = 1.0;
    options = simset('SrcWorkspace','current');
    sim('exsim',stop_time,options);
    analyze_result_func(bit_error,errdiff); 

    k = find(bit_error == 1, 1, 'last' )+1;
    if k <= 9000
        pass = true;
    else 
        pass = false;
    end
    cvrg = k;
    n_errors = sum(bit_error(bit_error ~= 0));

    errdiffx = errdiff(k:length(errdiff));
    MSE = mean(errdiffx.^2);
end