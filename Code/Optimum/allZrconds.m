function [ zinv,zinvs ] = allZrconds( chains,mode )
%ZRCONDS Summary of this function goes here
%   Detailed explanation goes here



zinv=NaN(size(chains));
zinvs=zinv;

for i=1:numel(chains)
    [zinv(i),zinvs(i)]=Zrconds(chains(i).s,chains(i).qv,mode);
end

end

