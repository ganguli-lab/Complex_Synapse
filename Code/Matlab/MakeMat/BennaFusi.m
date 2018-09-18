function [ Wp,Wm,w ] = BennaFusi( numvar,numlevel,m,varargin )
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

persistent p
if isempty(p)
    p=inputParser;
    p.FunctionName='BennaFusi';
    p.StructExpand=true;
    p.KeepUnmatched=false;
    p.addRequired('numvar',@(x)validateattributes(x,{'numeric'},{'scalar','integer','positive'},'BennaFusi','numvar',1))
    p.addRequired('numlevel',@(x)validateattributes(x,{'numeric'},{'scalar','even','positive'},'BennaFusi','numlevel',2))
    p.addRequired('m',@(x)validateattributes(x,{'numeric'},{'scalar','>',1},'BennaFusi','m',3))
    p.addOptional('decay',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'BennaFusi','decay',4));
    p.addOptional('rdt',1,@(x)validateattributes(x,{'numeric'},{'scalar','positive'},'BennaFusi','rdt',5));
end
p.parse(numvar,numlevel,m,varargin{:});
r=p.Results;
if any(strcmp('decay',p.UsingDefaults))
    r.decay=r.m^(-2*r.numvar+1);
end

% Each state can be described by a vector of levels, one element for each
% variable, in the range 0:r.numlevel-1. This can be converted to a state #:
% state # = sum_{i=1:numvar} level_i numlevel^(numvar-i)

numstates=r.numlevel^r.numvar;

Wp=zeros(numstates);
Wm=Wp;
w=ones(numstates/2,1);
w=[-w;w];

%state# = levelvec * levelsToState + 1;
%levelVec = BinVec(state# - 1, numlevel, numvar)
levelsToState = r.numlevel.^(r.numvar-1:-1:0)';


for i=1:numstates
    %state we jump from, row vector of levels of each variable
    fromst=wrev(BinVec(i-1,r.numlevel,r.numvar)); 
    du=diff(fromst);
    %state we jump to for potentiation, vector of levels of each variable
    %(real number, deal with integer part and probs later)
    tostp = fromst + (r.rdt/4)* r.m.^(-2*(1:r.numvar)+1) .* ( -r.m * [0 du] + [du 0]);
    %now deal with last var (decay to centre piece)
    tostp(end) = tostp(end) - r.rdt * r.decay * (fromst(end)-(r.numlevel-1)/2);
    %now for depression:
    tostm=tostp;
    %now deal with first var (pot vs dep piece):
    tostp(1) = tostp(1) + r.rdt;
    tostm(1) = tostm(1) - r.rdt;
    %prevent going over edge
    tostp = max(min(tostp,r.numlevel-1),0);
    tostm = max(min(tostm,r.numlevel-1),0);
    %now deal with integer part and probs
    for j=1:2^r.numvar
        %do we pick the level above or below tost(p/m)?
        whichjump=BinVec(j-1,2,r.numvar);
        %state we jump to
        jumpp=min(floor(tostp)+whichjump,r.numlevel-1);
        jumpm=min(floor(tostm)+whichjump,r.numlevel-1);
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

