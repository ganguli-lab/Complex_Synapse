function [ twochains ] = TwoModelDoubleEnv( chains,sc_ind,Ac )
%twochains=TWOMODELDOUBLEENV(chains,sc_ind,Ac) constrained envelope from
%population of two types of synapse
%   chains = struct from NumLaplaceBndChains

if Ac > chains(sc_ind).A
    error(['Ac = ' num2str(Ac) ' > env = ' num2str(chains(sc_ind).A)]);
end

n=length(chains(1).qv)/2;

srange=[chains.s];

twochains(1,length(chains)) = struct('s',[],'qv',[],'frac',[],'A',[],'snrb',[],'sc',srange(sc_ind),'Ac',Ac);

for i=1:length(chains)
    twochains(i).sc = srange(sc_ind);
    twochains(i).Ac = Ac;
    twochains(i).s = srange(i);
    twochains(i).qv = [chains(i).qv ; chains(sc_ind).qv];
    model=SynapseMemoryModel.Build(@MakeMultistate,0.5,{chains(i).qv(1:n),chains(i).qv(n+1:end)});
    Aother=model.SNRlaplace(srange(sc_ind));
    if Aother > Ac
        twochains(i).frac=[Ac/Aother ; 0];
        twochains(i).snrb=twochains(i).frac(1)*model.SNRrunAve(1./srange);
    else
        twochains(i).frac=[chains(sc_ind).A-Ac ; Ac-Aother]/(chains(sc_ind).A-Aother);
        twochains(i).snrb=twochains(i).frac(1)*model.SNRrunAve(1./srange);
        model=SynapseMemoryModel.Build(@MakeMultistate,0.5,{chains(sc_ind).qv(1:n),chains(sc_ind).qv(n+1:end)});
        twochains(i).snrb=twochains(i).snrb+twochains(i).frac(2)*model.SNRrunAve(1./srange);
    end
    twochains(i).A = twochains(i).snrb(i)/srange(i);
end



end

