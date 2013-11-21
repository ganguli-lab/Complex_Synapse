function [ newWp,newWm,neww,ix ] = SortByWtEtaS( Wp,Wm,w,fp,s )
%[newW,neww,ix]=SORTBYETA(M,w,fp) Put states in order of decreasing eta^+
%   M = transition rates, {Mpot,Mdep}
%   w = Weights of states (+/-1)
%   ix=sort order
%   newW=W(ix,ix)
%   neww=w(ix)

error(CheckSize(Wp,@ismat));%matrix
error(CheckSize(Wp,@issquare));%square
error(CheckSize(Wm,@samesize,'samesize(Wp)',Wp));%also square matrix of same size
error(CheckSize(fp,@isscalar));
error(CheckValue(fp,@inrange,'inrange(0,1)',0,1));%fp in [0,1]
error(CheckSize(w,@iscol));
error(CheckValue(w,@(x) all(x.^2==1),'all w = +/-1'));
error(CheckSize(w,@samelength,'samelength(Wp)',Wp));%same size

[neww,ix]=sort(w,'ascend');
newWp=Wp(ix,ix);
newWm=Wm(ix,ix);

Zinvs=s*eye(length(Wp)) + ones(length(Wp)) - fp*newWp - (1-fp)*newWm;
deta=-Zinvs\w;

wchange=[0 find(diff(w)) length(w)];

ix=[];
for i=1:length(wchange)-1
    [~,dix]=sort(deta(wchange(i)+1:wchange(i+1)),'descend');
    ix=[ix dix'+wchange(i)];
end

newWp=newWp(ix,ix);
newWm=newWm(ix,ix);
neww=neww(ix);


end