function [ newqv,A ] = ModelOptChainL( qv,sm,varargin)
%[newqv]=MODELOPTCHAINL(qv,tm) run gradient descent on model
%   sm = Laplace param value
%   qv = nearest neighbour transition rates

UseDerivs=false;
Algorithm='interior-point';
Display='off';
TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
DispExit=false;
varargin=assignApplicable(varargin);

qmin=zeros(size(qv));
qmax=ones(size(qv));

options = optimset('Algorithm',Algorithm,'Display',Display,...
    'TolFun', TolFun,...  % termination based on function value (of the derivative)
    'TolX', TolX,...
    'TolCon',TolCon,...
    'MaxIter',MaxIter, ...
    'largescale', 'on', ...
    varargin{:});

if UseDerivs
    options = optimset(options,'GradObj','on');
    [newqv,A,ef] = fmincon(@(y)-OptFunGradChainL(y,sm),qv,...
        [],[],[],[],qmin,qmax,[],... 
        options);
else
    [newqv,A,ef] = fmincon(@(y)-OptFunChainL(y,sm),qv,...
        [],[],[],[],qmin,qmax,[],... 
        options);
end
A=-A;

if DispExit
        switch ef
            case 1
                disp('Finished.  Magnitude of gradient smaller than the TolFun tolerance.');
%                 do_accept_opt_stop = true;
            case 2
                disp('Finished.  Change in x was smaller than the TolX tolerance.');
%                 do_accept_opt_stop = true;
            case 3 
                disp('Finished.  Change in the objective function value was less than the TolFun tolerance.');
%                 do_accept_opt_stop = true;
            case 5
                disp('Finished.  Predicted decrease in the objective function was less than the TolFun tolerance.');
%                 do_accept_opt_stop = true;
            case 0
                disp('Finished.  Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.');
%                 if do_topo_map
%                     do_accept_opt_stop = true;  % Hard to see how the fval will be less than fval_tol here, but anyways...
%                 else
%                     do_accept_opt_stop = false;
%                 end
            case -1
                disp('Finished.  Algorithm was terminated by the output function.');
                assert ( false, 'Still not sure what this case is.');
            otherwise
                assert ( false, 'New exit condition out of the fminunc optimizer front-end.');
        end
end

end

