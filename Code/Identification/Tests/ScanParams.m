
% tolfunvals=[10 5 1 0.5 0.1 0.05 0.01];
% tolxvals=[1e-2 5e-3 1e-3 5e-4 1e-4 5e-5 1e-5];
% crossvals=2:5;
% llincvals=[0 0.1 0.25 0.5 1 2 5 10];
NumSamples=[100 200 400 800];
priorCcoeff=[1 5 10 20 50];


scan_results(length(priorCcoeff),length(NumSamples))=...
    struct('num_events',[],'prob_st',[],'KL',[],'Ln',[]);

% options=SynapseOptimset('MaxStates',modelobj.NumStates/2+1,'ModelDiff','Ln');
% optimOptions=optimoptions('fmincon','Algorithm','interior-point','Display','off');
synapseOptions=SynapseOptimset('MaxStates',6);

for i=1:length(priorCcoeff)
    for j=1:length(NumSamples)
%         synapseOptions.TolFunChange=tolfunvals(i);
%         synapseOptions.TolX=tolxvals(j);
%         synapseOptions.MinLogLikeInc=llincvals(j);
        synapseOptions.PriorCcoeff=priorCcoeff(i);
        synapseOptions.NumSample=NumSamples(j);
        [scan_results(i,j).num_events,scan_results(i,j).prob_st,scan_results(i,j).KL,scan_results(i,j).Ln]=...
            TryParams( modelobj,synapseOptions,'crossval',1,'doDist',false );
    end
end

clear i j;