function [ Pertp,Pertm ] = ShortcutPert( W,start,len )
%[PERTP,PERTM] = SHORTCUTPERT(W,start,len) Adds shortcuts while preserving
%equilibrium probabilities
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   W = matrix we're permuting
%   start, = start of shortcut
%   len = length of shortcut

p=EqProb(W);
D=diag(1./p);

Pertp=zeros(size(W));

Pertp(start,start+len)=1;
for i=1:len
    Pertp(start+i-1,start+i)=-1;
end

Pertp=Stochastify(Pertp);
Pertm=Stochastify(Pertp');

Pertp=D*Pertp;
Pertm=D*Pertm;

end

