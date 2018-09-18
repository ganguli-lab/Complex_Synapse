function [ x ] = VORprettyfier( x0,builder_h,n,lambda,dS,varargin )
%VORPRETTYFIER Summary of this function goes here
%   Detailed explanation goes here

lbnds=zeros(size(x0));
ubnds=ones(size(x0));

lbnds(end-1:end)=[5 5];
ubnds(end-1:end)=[500 500];

options=optimset(varargin{:});

x=fmincon(@(y) VORprettyness(y,builder_h,n,lambda,dS), x0,...
    [],[],[],[],lbnds,ubnds,[],options);


end

