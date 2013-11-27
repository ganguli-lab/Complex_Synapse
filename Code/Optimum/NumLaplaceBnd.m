function [ chains ] = NumLaplaceBnd( srange,nstates,trange,varargin )
%chains=NUMLAPLACEBND(srange,nstates,trange) numeric laplace bound
%   chains  = struct array (size=[1 length(srange)])
%   srange  = values of Laplace parameter at which we maximise
%   nstates = number of states in chain
%   trange  = values of time for snr curve
%   chains.s   = value of Laplace parameter at which we optimised
%   chains.qv  = nearest neighbour transitions of optimal model
%   chains.A   = value of Laplace transform at s for optimal model
%   chains.snr = snr curve of optimal model

chains(1,length(srange))=struct('s',[],'qv',[],'A',[],'snr',[]);

for i=1:length(srange)
    chains(i).s=srange(i);
    [chains(i).qv,chains(i).A]=FindOptChainL(srange(i),nstates,50,varargin{:});
    [Wp,Wm,w]=MakeSMS(chains(i).qv);
    chains(i).snr=SNRcurve(trange,Wp,Wm,0.5,w);
    disp([int2str(i) '/' int2str(length(srange))]);
end


end

