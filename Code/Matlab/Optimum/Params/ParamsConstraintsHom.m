function [ A,b ] = ParamsConstraintsHom( n )
%[A,b]=PARAMSCONSTRAINTSHOM(n) Constrtaints on parameters from stochasticity
%   n = length(Wp)
%   Constaints: A.x <= b
%   From:           Wpm_ij >= 0, for all i ~= j
%        sum_{j~=i} Wpm_ij <= 1, for all i;
%                     Q_ij >= 0, for all i ~= j




A=zeros(2*n,3*n*(n-1));

for i=1:n
    A(i,(1+(i-1)*(n-1)):(i*(n-1)))=1;
    A(n+i,(n*(n-1)+1+(i-1)*(n-1)):(n*(n-1)+i*(n-1)))=1;
end

% A =[A; -eye(3*n*(n-1))];
% 
% b=[ones(2*n,1);zeros(3*n*(n-1),1)];
b=ones(2*n,1);



end

