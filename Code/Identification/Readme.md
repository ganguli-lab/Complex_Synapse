Synapse Identification
======================

This requires the Matlab File Exchange package mmx.
You will also need to build the mex files in folder BWmex

The main classes are SynapseIdModel, SynapsePlastSeq and SynapsePlastSeqSim

SynapseIdModel
--------------

This stores the complex synapse model that is the result of a fit 
or being used for a simulation.

*SynapseIdModel.M is a cell row of Markov matrices (discrete time) 
corresponding to each type of plasticty (e.g. potentiation and depression).
*SynapseIdModel.Initial is a row vector of probabilites for the initial state.
*SynapseIdModel.w is a column vector of synaptic weight labels for each state
SynapseIdModel.fp is a row vector of probabilities of each type of plasticty
being used, excluding the last one.

To create a ground truth model to use for simulation:
    modelobj=SynapseIdModel.Build(@SerialBuilder,{n,q},fp);
which creates a Serial synapse, with n states and transition probabilites q.

To create a random model as a starting guess for model fitting:
    guessmodel=SynapseIdModel.Rand(w);
which creates a model with length(w) states and random (normalised) Initial and M.
Optional parameter value pairs (with defaults): 
*'NumPlastTypes' number of typse of plasticity (2),
*'ScaleW' multiplies elelments (1),
*'sparsity' probability that each elemnt is non-zero (1).

SynapsePlastSeq
---------------

This stores the results of an experiment where a synapse is subjected to 
a sequence of different plasticity events and the synapstic weight is 
measured before each event.

*SynapsePlastSeq.potdep is row of plasticity types (index of SynapseIdModel.M)
used at each event.
*SynapsePlastSeq.readouts is a row of synaptic weight identifiers 
(1:number of values that SynapseIdModel.w can take) observed at each event.

To create an object for one sequence:
    seqobj=SynapsePlastSeq('potdep',potdep,'readouts',readouts);
where potdep and readouts are row vectors of the same length.

To store the results of many sequences, an array of SynapsePlastSeq can be 
used with different numbers of events in each one if needed.

SynapsePlastSeqSim
------------------

This is inherited from SynapsePlastSeq. It stores the results of a 
simulated experiment from a ground truth SynapseIdModel.

*SynapsePlastSeqSim.stateseq is a rowvector of synaptic state numbers 
before each plasticity event (same length as SynapsePlastSeq.potdep).

To simulate from SynapseIdModel modelobj:
    simobj=modelobj.Simulate(rand(2,t,m,n,...));
which creates an array of SynapsePlastSeqSim with size [m,n,...] and t events
in each one.


Fittng functions
----------------

Create an options struct with:
    options=SynapseOptimset('prameter',value,...);
read the comments of SynapseOptimset for details.

If you know the number of states needed:
    [fitmodel,fval,exitflag,output]=FitSynapse(simobj,guessmodel,options);
fitmodel is the result of the fit. fval is the log likelihood or log posterior.
If you only want to fit SynapseIdModel.M, use FitSynapseM.
If you only want to fit SynapseIdModel.Initial, use FitSynapseInit.

If you dont know the number of states needed:
    FitSynapseSize(simobj,options);
This will start with one state for each value of synaptic weight, 
then try adding one state for each synaptic weight and see if the 
likelihood/posterior increases.



