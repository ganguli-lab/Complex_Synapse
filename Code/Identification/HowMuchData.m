function [ num_data ] = HowMuchData( truemodel,varargin )
%num_data=HOWMUCHDATA(truemodel) How much data do we need to identify
%SynapseIdModel truemodel?
%   num_data = # plasticity events

DataPerChunk=50;
TestDataPerChunk=20;
KLdivThresh=1e-2;
fp=0.5;
Display=true;
StartLen=0;
MaxLen=512;
varargin=assignApplicable(varargin);

StartLen=floor(StartLen/2);

fitsim=truemodel.Simulate(fp,rand(2,DataPerChunk,StartLen));
testsim=truemodel.Simulate(fp,rand(2,TestDataPerChunk,StartLen));


while length(fitsim)<MaxLen
    DoubleData;
    fitmodel=FitSynapseSize(fitsim,testsim,'Display',false,varargin{:});
    if truemodel.SameSizes(fitmodel)
        [~,kl]=truemodel.KLdivs(fitmodel);
        if mean(kl)/truemodel.NumStates < KLdivThresh;
            break;
        end%if KLdivThresh
    end%if SameSizes
end

num_data = fitsim.NumT + testsim.NumT;

    function DoubleData
        extrasim=truemodel.Simulate(fp,rand(2,DataPerChunk,max(1,length(fitsim))));
        fitsim=[fitsim extrasim];
        %
        extrasim=truemodel.Simulate(fp,rand(2,TestDataPerChunk,max(1,length(testsim))));
        testsim=[testsim extrasim];
        %
        if Display
            disp(length(fitsim));
        end%if Display
    end%function DoubleData

end

