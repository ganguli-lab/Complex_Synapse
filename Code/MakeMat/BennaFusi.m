function [ Wp,Wm,w ] = BennaFusi( numvar,numlevel,m,decay,rdt )
%[Wp,Wm,w]=BENNAFUSI(numvar,numlevel,m,rdt) Benna-Fusi variable diffusion
%model.
%   WP = potentiation transition rates
%   WM = depression transition rates
%   w  = Weights of states (+/-1)
%   numvar   = number of varibales (min=2)
%   numlevel = number of levels that each variable can take (even)
%   m        = sqrt of ratio of adjacent timescales
%   decay    = decay rate of last state (default = m^(-2*numvar+1))
%   rdt      = rate X time-step (default=1)

error(CheckSize(numvar,@isscalar));
error(CheckValue(numvar,@isint));
error(CheckSize(numlevel,@isscalar));
error(CheckValue(numlevel,@(x) isint(x/2),'even'));
error(CheckSize(m,@isscalar));
error(CheckValue(m,@(x) x>1,'> 1'));
existsAndDefault('decay', m^(-2*numvar+1));
error(CheckSize(decay,@isscalar));
existsAndDefault('rdt',1);
error(CheckSize(rdt,@isscalar));

% Each state can be described by a vector of levels, one element for each
% variable, in the range 0:numlevel-1. This can be converted to a state #:
% state # = sum_{i=1:numvar} level_i numlevel^(numvar-i)

numstates=numlevel^numvar;

Wp=zeros(numstates);
Wm=Wp;
w=ones(numstates/2,1);
w=[-w;w];

%state# = levelvec * levelsToState + 1;
%levelVec = BinVec(state# - 1, numlevel, numvar)
levelsToState = numlevel.^(numvar-1:-1:0)';


for i=1:numstates
    %state we jump from, row vector of levels of each variable
    fromst=wrev(BinVec(i-1,numlevel,numvar)); 
    du=diff(fromst);
    %state we jump to for potentiation, vector of levels of each variable
    %(real number, deal with integer part and probs later)
    tostp = fromst + (rdt/4)* m.^(-2*(1:numvar)+1) .* ( -m * [0 du] + [du 0]);
    %now deal with last var (decay to centre piece)
    tostp(end) = tostp(end) - rdt * decay * (fromst(end)-(numlevel-1)/2);
    %now for depression:
    tostm=tostp;
    %now deal with first var (pot vs dep piece):
    tostp(1) = tostp(1) + rdt;
    tostm(1) = tostm(1) - rdt;
    %prevent going over edge
    tostp = max(min(tostp,numlevel-1),0);
    tostm = max(min(tostm,numlevel-1),0);
    %now deal with integer part and probs
    for j=1:2^numvar
        %do we pick the level above or below tost(p/m)?
        whichjump=BinVec(j-1,2,numvar);
        %state we jump to
        jumpp=min(floor(tostp)+whichjump,numlevel-1);
        jumpm=min(floor(tostm)+whichjump,numlevel-1);
        %prob of making that jump
        probp=1-abs(tostp-jumpp);
        probm=1-abs(tostm-jumpm);
        %now put it in the Markov matrix
        Wp(i,jumpp*levelsToState+1)=prod(probp);
        Wm(i,jumpm*levelsToState+1)=prod(probm);
    end%for j
end%for i

%now deal with diagonal elements
Wp=StochastifyC(Wp);
Wm=StochastifyC(Wm);


end

