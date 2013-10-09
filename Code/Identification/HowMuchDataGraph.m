function [ data_s ] = HowMuchDataGraph( builder,param,varargin )
%data_s=HOWMUCHDATAGRAPH Summary of this function goes here
%   data_s = struct(n,data)
%   builder = function handle [Wp,Wm,w]=builder(n,param)


MinN=2;
MaxN=12;
StartLen=2;
varargin=assignApplicable(varargin);


data_s = struct('n',MinN:2:MaxN,'data',NaN(1,(MaxN-MinN)/2+1));

for i=0:(MaxN-MinN)/2
    n=MinN+2*i;
    modelobj=SynapseIdModel.Build(builder,0.5,n,param);
    try
        data_s.data(i+1)=HowMuchData(modelobj,'MaxStates',n/2+1,'StartLen',StartLen,varargin{:});
    catch ME
        disp(['n=' int2srt(n) '. ' ME.message]);
    end
    if ~isnan(data_s.data(i+1))
        StartLen=data_s.data(i+1);
    end
end

end

