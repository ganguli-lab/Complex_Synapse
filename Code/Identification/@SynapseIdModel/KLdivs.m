function [ initdiv,Mdivs ] = KLdivs( obj1,obj2 )
%[initdiv,Mdivs]=KLDIVS(obj1,obj2) Kullback-Leibler divergences from obj1
%to obj2
%   initdiv  = KL div of obj.Initial
%   Mdivs(i) = KL div of obj.M{i}

initdiv=KLdiv(obj1.Initial,obj2.Initial);
Mdivs=cellfun(@KLdiv,obj1.M,obj2.M);
Mdivs=Mdivs/obj1.NumStates;

end

