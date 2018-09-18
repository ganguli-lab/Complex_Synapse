function [ A,b ] = HomChainConstr( numstates, fp )
%[A,b]=HOMCHAINCONSTR(numstates,fp) contraints on transition probs for chain
%with homeostatic plasticity
%   constraints are A * q <= b
%   acts on q = [q_pot_inc q_pot_dec q_dep_inc q_dep_dec]'
fm=1-fp;

zermat=zeros(numstates-1);
onemat=eye(numstates-1);

zervec=zeros(numstates-1,1);
onevec=ones(numstates-1,1);

%positivity
A=-eye(4*(numstates-1));
b=zeros(4*(numstates-1),1);

%fp*Mp(i,i+1)-fm*Mm(i,i+1) >= 0
A =[A; [-fp*onemat zermat fm*onemat zermat]];
b=[b; zervec];

%fp*Mp(i,i-1)-fm*Mm(i,i-1) <= 0
A =[A; [zermat fp*onemat zermat -fm*onemat]];
b=[b; zervec];

%fp*Mp(i,i+1)-fm*Mm(i,i+1) <= fp
A =[A; [fp*onemat zermat -fm*onemat zermat]];
b=[b; fp*onevec];

%fp*Mp(i,i+1)-fm*Mm(i,i+1) >= -fm
A =[A; [zermat -fp*onemat zermat fm*onemat]];
b=[b; fm*onevec];



end

