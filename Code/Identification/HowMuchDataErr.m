function [ err_data ] = HowMuchDataErr( truemodel,varargin )
%err_data=HowMuchDataErr(truemodel) How much data do we need to identify
%SynapseIdModel truemodel?
%   err_data = struct(num_data,KLdiv_mean,KLdiv_err,L2_mean,L2_err)

DataPerChunk=200;
fp=0.5;
Display=true;
StartLen=1;
MaxLen=512;
NumReps=10;
SuccessFrac=0.5;
varargin=assignApplicable(varargin);
TestDataPerChunk=DataPerChunk;
varargin=assignApplicable(varargin);



err_data = struct('num_data',zeros(1,ceil(log2(MaxLen/StartLen))));
err_data.KLdiv_mean=NaN(size(err_data.num_data));
err_data.KLdiv_err=err_data.KLdiv_mean;
err_data.L2_mean=err_data.KLdiv_mean;
err_data.L2_err=err_data.KLdiv_mean;

StartLen=StartLen/2;

for i=1:ceil(log2(MaxLen/StartLen))
    StartLen=StartLen*2;
    if Display
        disp(StartLen);
    end
    klv=NaN(1,NumReps);
    L2v=klv;
    for j=1:NumReps
        fitsim=truemodel.Simulate(fp,rand(2,DataPerChunk,StartLen));
        testsim=truemodel.Simulate(fp,rand(2,TestDataPerChunk,StartLen));
        %
        err_data.num_data(i) = fitsim.NumT + testsim.NumT;
        %
        fitmodel=FitSynapseSize(fitsim,testsim,'Display',false,varargin{:});
        if truemodel.SameSizes(fitmodel)
            [~,kl]=truemodel.KLdivs(fitmodel);
            [~,L2]=truemodel.LnNorm(2,fitmodel);
            klv(j)=mean(kl)/truemodel.NumStates;
            L2v(j)=sum(L2);
        end%if SameSizes
    end%for j
    klv(isnan(klv))=[];
    L2v(isnan(L2v))=[];
    if length(klv)>=SuccessFrac*NumReps
        err_data.KLdiv_mean(i)=mean(klv);
        err_data.KLdiv_err(i)=std(klv)/sqrt(length(klv));
        err_data.L2_mean(i)=mean(L2v);
        err_data.L2_err(i)=std(L2v)/sqrt(length(L2v));
    end%if length(klv
end


end

