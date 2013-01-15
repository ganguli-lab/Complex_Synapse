function MakeVORcomparison( paramWT,paramKn,df,n,modelname )
%MAKEVORCOMPARISON(paramp,paramm,df,n,modelname) comparison of VOR learning
%curves for a particular synaptic model
%   paramWT   = parameter used for wild-type (pot&dep) and knockout (pot)
%   paramKn   = parameter used for knockout (dep)
%   df        = change in f+ vor gain increase learning
%   n         = number of internal synaptic states
%   modelname = name of synaptic model (multistate or cascade)

if strcmpi(modelname,'multistate')
    [Wp,Wm,w]=MakeSMS(paramWT*ones(1,n-1));
    [~,WmKn]=MakeSMS(paramKn*ones(1,n-1));
elseif strcmpi(modelname,'cascade')
    [Wp,Wm,w]=CascadeOriginal(paramWT,paramWT,n);
    [~,WmKn]=CascadeOriginal(paramWT,paramKn,n);
else
    error(['invalid model name: ' modelname]);
end


VORcomparison(Wp,Wm,w,0:1:80,df,40,Wp,WmKn,'LineWidth',2);

end

