function [ Wp,Wm,w ] = MakeSMSsq( qvec,fp )
%[WP,WM]=MAKESMS(QVEC) Transition rates for SMS toplogy
%   QVEC=adjacent transiton rates for potentiation, row vector

assert(isrow(qvec));

[ Wp,Wm,w ] = MakeSMS( qvec );

W=fp*Wp+(1-fp)*Wm;
W=-W*W;
[Wp,Wm]=TriangleDcmp(W,fp);



end

