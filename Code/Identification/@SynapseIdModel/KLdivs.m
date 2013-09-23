function [ initdiv,Mdivs ] = KLdivs( obj1,obj2 )
%[initdiv,Mdivs]=KLDIVS(obj1,obj2) Kullback-Leibler divergences from obj1
%to obj2
%   initdiv  = KL div of obj.initial
%   Mdivs(i) = KL div of obj.M{i}

initdiv=KLdiv(obj1.initial,obj2.initial);
Mdivs=cellfun(@KLdiv,obj1.M,obj2.M);

end

