function [ initdiv,Mdivs ] = KLdivs( obj1,obj2 )
%[initdiv,Mdivs]=KLDIVS(obj1,obj2) Kullback-Leibler divergences from obj1
%to obj2
%   initdiv  = KL div of obj.Initial
%   Mdivs(i) = KL div of obj.M{i}

initdiv=KLcalc(obj1.Initial,obj2.Initial);
Mdivs=cellfun(@KLcalc,obj1.M,obj2.M);
Mdivs=Mdivs/obj1.NumStates;


    function div=KLcalc(A,B)
        B=B(A>0);
        A=A(A>0);

        div=sum(A.*log(A./B));
    end

end

