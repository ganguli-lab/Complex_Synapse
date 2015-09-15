function [ newenv ] = MixDoubleEnv( env )
%MIXDOUBLEENV Summary of this function goes here
%   Detailed explanation goes here

snrb=reshape([env.mats.snrb],[],length(env.mats));

[SNRbenv,ix]=max(snrb,[],2);


newenv=env;
newenv.mats=env.mats(ix);
newenv.SNRbenv=SNRbenv;
SNRbenv=num2cell(SNRbenv);
[newenv.mats.A]=deal(SNRbenv{:});


end

