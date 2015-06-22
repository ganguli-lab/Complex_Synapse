function [ msg ] = ExitFlagMsg( ef )
%msg=EXITFLAGMSG(ef) Convert exit flag from FMINCON to messgae
%   Detailed explanation goes here

        switch ef
            case 1
                msg=('Finished.  Magnitude of gradient smaller than the TolFun tolerance.');
%                 do_accept_opt_stop = true;
            case 2
                msg=('Finished.  Change in x was smaller than the TolX tolerance.');
%                 do_accept_opt_stop = true;
            case 3 
                msg=('Finished.  Change in the objective function value was less than the TolFun tolerance.');
%                 do_accept_opt_stop = true;
            case 5
                msg=('Finished.  Predicted decrease in the objective function was less than the TolFun tolerance.');
%                 do_accept_opt_stop = true;
            case 0
                msg=('Finished.  Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.');
%                 if do_topo_map
%                     do_accept_opt_stop = true;  % Hard to see how the fval will be less than fval_tol here, but anyways...
%                 else
%                     do_accept_opt_stop = false;
%                 end
            case -1
                msg=('Finished.  Algorithm was terminated by the output function.');
            otherwise
                assert ( false, 'New exit condition out of the fminunc optimizer front-end.');
        end

end

