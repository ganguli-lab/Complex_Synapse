function printVORfigs( prefix,paramWT,paramKO,df,T_train,T_pre,varargin )
%PRINTVORFIGS Summary of this function goes here
%   Detailed explanation goes here

paramPot=paramWT;
pooled=false;
if strncmpi(prefix,'cascade',3)
    builder_h=@CascadeBuilder;
    n=10;
elseif strncmpi(prefix,'nonuni',3)
    builder_h=@NonuniBuilder;
    n=10;
elseif strncmpi(prefix,'serial',3)
    builder_h=@SerialBuilder;
    n=10;
elseif strncmpi(prefix,'multistate',3)
    builder_h=@MultistateBuilder;
    n=10;
elseif strncmpi(prefix,'binary',3)
    builder_h=@SerialBuilder;
    n=2;
elseif strncmpi(prefix,'pooled',3)
    builder_h=@PooledBuilder;
%     pooled=true;
    paramPot=paramWT(3)*[1 1];
    paramWT(3)=[];
    n=7;
end



vexpt=VORbuilder(builder_h,n,paramPot,paramWT,paramKO,0.5,0.5+df,0.5-df,T_train,T_pre);
vexpt.PrintFigs(prefix);


end

