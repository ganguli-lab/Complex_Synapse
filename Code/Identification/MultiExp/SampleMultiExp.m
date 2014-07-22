function [ samples ] = SampleMultiExp( params , num_samples)
%SAMPLEMULTIEXP Summary of this function goes here
%   Detailed explanation goes here


n=ceil(length(params)/2);
tau=exp(-params(1:n));
c=params(n+1:end);
c=[c;1-sum(c)];

ccum=cumsum(c);

samples=zeros(1,num_samples);

for i=1:num_samples
    testnum=rand;
    whichexp=find(testnum<ccum,1,'first');
    pd=makedist('Exponential','mu',tau(whichexp));
    samples(i)=random(pd,1);
end


end

