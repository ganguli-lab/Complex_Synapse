function [ initnrm,Mnrm ] = LnNorm( obj1,n,obj2 )
%[initnrm,Mnrm]=LNNORM(obj1,n,obj2) L_n norm of obj1 or obj1-obj2
%   initnrm  = L_n norm of obj.Initial
%   Mnrm(i) = L_n norm of obj.M{i}

existsAndDefault('n',2);

if exist('obj2','var')
    [ initnrm,Mnrm ] = LnNorm( obj1-obj2,n );
else
    if n==0
        initnrm=sum(obj1.Initial~=0);
        Mnrm=cellfun(@(x) sum(reshape(x~=0,1,[])),obj1.M);
    elseif isinf(n)
        initnrm=max(obj1.Initial);
        Mnrm=cellfun(@(x) max(reshape(x,1,[])),obj1.M);
    else
        initnrm=sum(abs(obj1.Initial).^n)^(1/n);
        Mnrm=cellfun(@(x) sum(reshape(abs(x),1,[]).^n)^(1/n),obj1.M);
    end
end

end

