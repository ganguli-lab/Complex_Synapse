function [ initnrm,Mnrm ] = LnNorm( obj1,obj2 )
%[initnrm,Mnrm]=LNNORM(obj1,obj2) L_n norm of obj1 or obj1-obj2
%   initnrm  = L^n norm of obj.Initial
%   Mnrm(i) = combination of L^n norm of rows of obj.M{i}
%   n = obj1.NormPower
%   obj1.NormRows used to combine rows (e.g. @max, @mean, @norm)


if exist('obj2','var')
    [ initnrm,Mnrm ] = LnNorm( obj1-obj2 );
else
    if obj1.NormPower==0
        initnrm=sum(obj1.Initial~=0);
        Mnrm=cellfun(@(x) mean(sum(x~=0,2)),obj1.M);
        Mnrm=Mnrm/obj1.NumStates;
    elseif isinf(obj1.NormPower)
        initnrm=max(obj1.Initial);
        Mnrm=cellfun(@(x) mean(max(x,[]),2),obj1.M);
    else
        initnrm=sum(abs(obj1.Initial).^obj1.NormPower)^(1/obj1.NormPower);
        Mnrm=cellfun(@(x) obj1.NormRows(sum(abs(x).^obj1.NormPower,2).^(1/obj1.NormPower)),obj1.M);
    end
end

end

