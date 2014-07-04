
tolfunvals=[10 5 1 0.5 0.1 0.05 0.01];
tolxvals=[1e-2 5e-3 1e-3 5e-4 1e-4 5e-5 1e-5];

scan_results(length(tolfunvals),length(tolxvals))=...
    struct('num_events',[],'prob_st',[],'KL',[],'Ln',[]);

options=SynapseOptimset('MaxStates',modelobj.NumStates/2,'ModelDiff','Ln');

for i=1:length(tolfunvals)
    for j=1:length(tolxvals)
        options.TolFunChange=tolfunvals(i);
        options.TolX=tolxvals(j);
        [scan_results(i,j).num_events,scan_results(i,j).prob_st,scan_results(i,j).KL,scan_results(i,j).Ln]=...
            TryParams( modelobj,options,'doSize',false );
    end
end

clear i j;