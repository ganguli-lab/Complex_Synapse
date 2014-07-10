function [ params,negloglike,exitflag,output,lambda,grad,hessian ] = FitMultiExp( num_exp,data,optimOptions )
%FITMULTIEXP Summary of this function goes here
%   Detailed explanation goes here


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

