function [ params,negloglike,exitflag,output,lambda,grad,hessian ] = FitMultiExp( num_exp,data,optimOptions )
%[params,negloglike,exitflag,output,lambda,grad,hessian]=FITMULTIEXP(num_exp,data,optimOptions)
%max-likelihood fit for multi-exponential distribution of data
%   num_exp       = number of exponentials used for fit
%   data          = row of dwell times
%   optimOptions  = optimoptions for fmincon
%   params = column [q; c(1:end-1)]
%       c(end)=1-sum(c(1:end-1))
%       like = prod_data sum_{c,q} c*q*exp(-q*data)
%   negloglike = -log(like)
%exitflag,output,lambda,grad,hessian from fmincon

A=-eye(num_exp,2*num_exp-1);
b=zeros(num_exp,1);
guess=rand(2*num_exp-1,1);

if strcmp(optimOptions.GradObj,'on')
    if strcmp(optimOptions.Hessian,'user-supplied')
        fun=@MultiExpFunGradHess;
    else
        fun=@MultiExpFunGrad;
    end
else
    fun=@MultiExpFun;
end

[ params,negloglike,exitflag,output,lambda,grad,hessian ] = fmincon(@(x)fun(x,data),guess,A,b,[],[],[],[],[],optimOptions);
end

