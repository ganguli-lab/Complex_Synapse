function [ tf,ix,fail ] = CheckSplit( obj,varargin )
%[tf,ix]=obj.CHECKSPLIT check if SynapseMemoryModel, obj , consists of
%disconnected models
%   checks smallest singular value

[U,S,~]=svd(obj.GetZinv);

[smin,ix]=min(diag(S));
U=U(:,ix);
ix=U>0;

fail = smin < obj.SingValThresh && (all(ix) || all(~ix));

tf = smin < obj.SingValThresh && ~all(ix) && ~all(~ix);


end

