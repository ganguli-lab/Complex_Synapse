srange=10.^(-4:0.1:1);
trange=10.^(-1:0.1:3);


chains=NumLaplaceBndChain(srange,12,trange,'other','DispReps',true);
AenvChains=[chains.A];
Q=reshape([chains.qv],[],51);

save laplacebndChainA12 chains AenvChains Q srange trange;

clear chains AenvChains Q;

mats=NumLaplaceBnd(srange,12,trange,'other','DispReps',true);
AenvAll=[mats.A];

save laplacebnd12 mats AenvAll srange trange;