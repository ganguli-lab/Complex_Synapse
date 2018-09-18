function [ h ] = ShortcutMix( nstates,transprob,fp, start,len, pertvals,  varargin )
%H=PERTPLOT(WP,WM,PERTP,PERTM,FP,w,PERTVALS,CALCHANDLE,...) plot of area under SNR
%curve for a perturbation
%   H=plot handle
%   WP,WM = base transition rates for potentiation,depression
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   FP = fraction of potentiation events
%   w = Weights of states (+/-1)
%   PERTVALS = vector/matrix of perturbation multipliers for plot
%   CALCHANDLE = function handle, acts on (Wp,Wm,fp,w)
%   ... = passed to plot function

q=ones(1,nstates-1)*transprob;
[wp,wm]=MakeSMS(q);
W=fp*wp+(1-fp)*wm;
tau=MixTime(W);


[pertp,pertm]=ShortcutPert(W,start,len);
Pert=fp*pertp+(1-fp)*pertm;
        
% vals=zeros(3,length(pertvals));
vals=zeros(1,length(pertvals));

for i=1:numel(pertvals)
    R=(nstates*pertvals(i))/(transprob);
    vals(1,i)= tau/MixTime(W+pertvals(i)*Pert);
%     vals(2,i)=1/(1+R);
%     vals(3,i)=1+(len-1)*R;
end

h=plot(pertvals,vals,varargin{:});

end
