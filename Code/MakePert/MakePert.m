function [ Pertp,Pertm ] = MakePert( siz,row,col )
%[PERTP,PERTM] = MAKEPERT(SIZ,ROW,COL) Summary of this function goes here
%   PERTP,PERTM = Perturbation of transition rates for potentiation,depression
%   SIZ = size of matrices
%   row,col = row,col of W^+ we're perturbing

        Pertp=zeros(siz);
        Pertp(row,col)=1;
        Pertp(row,row)=-1;
        Pertm=rot90(Pertp,2);



end

