% EXAMPLE  An example which calls ucordic.m
%
% << Syntax >>
%    max_abs_error = example(p1,n1,p2,n2,Niter)

function max_abs_error = example(p1,n1,p2,n2,Niter)

    % Testing angles in the range from 2 to 5 degrees, 0.5 deg step
    
    run_init_step = false;
    
    angles_deg = (2:0.5:5);
    angles_rad = pi * angles_deg / 180; % Note: angles must be in radians in ucordic
    Nangles = length(angles_rad);

    fixptpar.xy_p = p1;
    fixptpar.xy_n = n1;
    fixptpar.z_p = p2;
    fixptpar.z_n = n2;    
    gaincomp = true; 
    
    override = false;
    verbose = false;
    
    max_abs_error = 0; 
    for n = 1:Nangles        
        % Tested angle. Mapping it to the range [-pi,pi].
        phi = mod(angles_rad(n)+pi,2*pi)-pi;        

        % Setting run_init_step = true above allows using any angle  
        % Works as explained in Section 2.5. of Intro III (1/2)
        if run_init_step
            if phi < 0
                d_init = -1;
            else
                d_init = +1;
            end
            x0 = 0;
            y0 = d_init;
            z0 = phi - d_init * pi / 2;
        else
            x0 = 1;
            y0 = 0;
            z0 = phi;
        end
	
        % Double precision "ground truth"
        cosval = cos(phi); 
        sinval = sin(phi);
        
        % Run CORDIC        
        [xN,yN,~,~] = ucordic('circular','rotation',x0,y0,z0,Niter,gaincomp,override,fixptpar,verbose);        
        fprintf('Angle: %.2f deg | Cos error | = %.6f | Sin error | = %.6f\n',angles_deg(n),abs(cosval-xN),abs(sinval-yN));
        
        % Update maximum observed absolute error
        emax = max(abs(cosval-xN),abs(sinval-yN));
        if emax > max_abs_error
            max_abs_error = emax;
        end
        
    end
    
end

