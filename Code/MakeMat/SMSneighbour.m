function [ Wp,Wm,w,Wpu,Wmu ] = SMSneighbour( len,fp,direct,shortcut,sparsity )
%SMSNEIGHBOUR Summary of this function goes here
%   Detailed explanation goes here


q=ones(1,len-1)*direct;
[Wpu,Wmu]=MakeSMS(q);
w=[-ones(len/2,1);ones(len/2,1)];
[Wp,Wm,w]=TriangleDcmp(fp*Wpu+(1-fp)*Wmu+RandTrans(len,sparsity)*shortcut,fp,w);

end

