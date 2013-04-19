function [ A,b ] = ParamsConstraints( n )
%[A,b]=PARAMSCONSTRAINTS(n) Constrtaints on parameters from stochasticity
%   n = length(Wp)
%   Constaints: A.x <= b
%   From:           Wpm_ij >= 0, for all i ~= j
%        sum_{j~=i} Wpm_ij <= 1, for all i;




A=zeros(2*n,2*n*(n-1));

for i=1:n
    A(i,(1+(i-1)*(n-1)):(i*(n-1)))=1;
    A(n+i,(n*(n-1)+1+(i-1)*(n-1)):(n*(n-1)+i*(n-1)))=1;
end

A =[A; -eye(2*n*(n-1))];

b=[ones(2*n,1);zeros(2*n*(n-1),1)];



end

