function [ env ] = TwoModelDoubleEnv( chains,sc_ind,Ac )
%env=TWOMODELDOUBLEENV(chains,sc_ind,Ac) constrained envelope from
%population of two types of synapse
%   chains = struct from NumLaplaceBndChains

if Ac > chains(sc_ind).A
    error(['Ac = ' num2str(Ac) ' > env = ' num2str(chains(sc_ind).A)]);
end


tau=1./[chains.s];
SNRbc=Ac*chains(sc_ind).s;

twochains(1,length(chains)) = struct('s',[],'qv',[],'frac',[],'A',[],'snrb',[],'sinds',[]);

for i=1:length(chains)
    twochains(i).s = chains(i).s;
    twochains(i).A = 0;
    for j=1:length(chains)
        for k=1:length(chains)
            tempchains=twochains(i);
            tempchains.qv = [chains(j).qv ; chains(k).qv];
            tempchains.sinds=[j k];

            Sij=chains(j).snrb(i);
            Sik=chains(k).snrb(i);
            Scj=chains(j).snrb(sc_ind);
            Sck=chains(k).snrb(sc_ind);
            
            if Scj > SNRbc && Sck > SNRbc
                x=Sij>Sik;
            elseif Scj < SNRbc && Sck < SNRbc
                x=NaN;
            elseif Scj > Sck && Sij > Sik
                x=1;
            elseif Scj < Sck && Sij < Sik
                x=0;
            else
                x=(SNRbc-Sck)/(Scj-Sck);
            end
            tempchains.frac = [x 1-x];
            tempchains.snrb = tempchains.frac * [chains(j).snrb; chains(k).snrb];
            tempchains.A = tempchains.snrb(i)/tempchains.s;
            
            if tempchains.A > twochains(i).A
                twochains(i)=tempchains;
            end
        end%for k
    end%for j
end

Aenv=[twochains.A];
Aenv=Aenv./tau;


env=struct('chains',twochains,'sc',chains(sc_ind).s,'Ac',Ac,'tau',tau,'SNRbenv',Aenv);

end

