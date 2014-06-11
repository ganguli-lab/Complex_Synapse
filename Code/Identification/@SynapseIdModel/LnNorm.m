function [ initnrm,Mnrm ] = LnNorm( obj1,n,obj2 )
%[initnrm,Mnrm]=LNNORM(obj1,n,obj2) L_n norm of obj1 or obj1-obj2
%   initnrm  = L_n norm of obj.Initial
%   Mnrm(i) = mean of L^n norm of rows of obj.M{i}

existsAndDefault('n',2);

if exist('obj2','var')
    [ initnrm,Mnrm ] = LnNorm( obj1-obj2,n );
else
    if n==0
        initnrm=sum(obj1.Initial~=0);
        Mnrm=cellfun(@(x) mean(sum(x~=0,2)),obj1.M);
        Mnrm=Mnrm/obj1.NumStates;
    elseif isinf(n)
        initnrm=max(obj1.Initial);
        Mnrm=cellfun(@(x) mean(max(x,[]),2),obj1.M);
    else
        initnrm=sum(abs(obj1.Initial).^n)^(1/n);
        Mnrm=cellfun(@(x) mean(sum(abs(x).^n,2).^(1/n)),obj1.M);
    end
end

end

