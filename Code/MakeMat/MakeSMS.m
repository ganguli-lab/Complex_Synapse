function [ Wp,Wm,w ] = MakeSMS( qvec )
%[WP,WM,w]=MAKESMS(QVEC) Transition rates for SMS toplogy
%   QVEC=adjacent transiton rates for potentiation, row vector

error(CheckSize(qvec,@isrow));

Wp=diag(qvec,1)-diag([qvec 0]);

if nargout>1
    Wm=rot90(Wp,2);
end

if nargout>2
    w=BinaryWeights(length(qvec)+1);
end

end

