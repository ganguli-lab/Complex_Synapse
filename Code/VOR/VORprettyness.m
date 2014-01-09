function [ f ] = VORprettyness( x,builder_h,n,lambda,dS )
%VORPRETTYNESS Summary of this function goes here
%   Detailed explanation goes here

x=num2cell(x);
vexpt=VORbuilder(builder_h,n,x{:});
S=vexpt.LearnSdata;

f=(S.WTnopre(end)-S.KOwithpre(end))^2 ...
    +(S.KOnopre(end)-S.WTwithpre(end))^2 ...
    +lambda*(S.WTnopre(end)-S.KOnopre(end)-dS*abs(S.WTnopre(end)-S.KOwithpre(end)))^2;

end

