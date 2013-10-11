function [ err_data ] = HowMuchDataErr( truemodel,varargin )
%err_data=HowMuchDataErr(truemodel) How much data do we need to identify
%SynapseIdModel truemodel?
%   err_data = struct(num_data,KLdiv,L2err)

DataPerChunk=200;
TestDataPerChunk=50;
fp=0.5;
Display=true;
StartLen=1;
MaxLen=512;
varargin=assignApplicable(varargin);



err_data = struct('num_data',zeros(1,ceil(log2(MaxLen/StartLen))));
err_data.KLdiv=NaN(size(err_data.num_data));
err_data.L2err=NaN(size(err_data.num_data));

StartLen=floor(StartLen/2);
i=0;
fitsim=truemodel.Simulate(fp,rand(2,DataPerChunk,StartLen));
testsim=truemodel.Simulate(fp,rand(2,TestDataPerChunk,StartLen));

while length(fitsim)<MaxLen
    DoubleData;
    i=i+1;
    fitmodel=FitSynapseSize(fitsim,testsim,'Display',false,varargin{:});
    err_data.num_data(i) = fitsim.NumT + testsim.NumT;
    if truemodel.SameSizes(fitmodel)
        [~,kl]=truemodel.KLdivs(fitmodel);
        [~,L2]=truemodel.LnNorm(2,fitmodel);
        err_data.KLdiv(i)=mean(kl)/truemodel.NumStates;
        err_data.L2err(i)=sum(L2);
    end%if SameSizes
end


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

