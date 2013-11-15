function [ chains ] = NumLaplaceBnd( srange,nstates,trange,varargin )
%NUMLAPLACEBND Summary of this function goes here
%   Detailed explanation goes here

chains(1,length(srange))=struct('s',[],'qv',[],'A',[],'snr',[]);

for i=1:length(srange)
    chains(i).s=srange(i);
    [chains(i).qv,chains(i).A]=FindOptChainL(srange(i),nstates,50,varargin{:});
    [Wp,Wm,w]=MakeSMS(chains(i).qv);
    chains(i).snr=SNRcurve(trange,Wp,Wm,0.5,w);
    disp([int2str(i) '/' int2str(length(srange))]);
end


end

