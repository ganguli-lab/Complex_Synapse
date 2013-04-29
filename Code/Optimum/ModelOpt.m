function [ newWp,newWm ] = ModelOpt( Wp,Wm,tm,varargin)
%[newWp,newWm]=MODELOPT(Wp,Wm,tm) run gradient descent on model
%   T = time value
%   WP = potentiation transition rates
%   WM = depression transition rates

UseDerivs=false;
Algorithm='interior-point';
Display='off';
TolFun=1e-6;
TolX=1e-10;
TolCon=1e-6;
MaxIter=1000;
varargin=assignApplicable(varargin);

n=length(Wp);
w=BinaryWeights(n);

[A,b]=ParamsConstraints(n);

x0 = Mats2Params(Wp,Wm);            % Starting guess 
options = optimset('Algorithm',Algorithm,'Display',Display,...
    'TolFun', TolFun,...  % termination based on function value (of the derivative)
    'TolX', TolX,...
    'TolCon',TolCon,...
    'MaxIter',MaxIter, ...
    'largescale', 'on', ...
    varargin{:});

if UseDerivs
    options = optimset(options,'GradObj','on');
    [x,~,ef] = fmincon(@(y)OptFunGrad(y,tm,0.5,w),x0,A,b,...
         [],[],[],[],[],... 
       options);
else
    [x,~,ef] = fmincon(@(y)OptFun(y,tm,0.5,w),x0,A,b,...
         [],[],[],[],[],... 
       options);
end
[Wp,Wm]=Params2Mats(x);

[~,~,ix]=SortByEta(0.5*Wp+0.5*Wm,w);
newWp=Wp(ix,ix);
newWm=Wm(ix,ix);

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

