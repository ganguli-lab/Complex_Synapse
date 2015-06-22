function [ tm ] = TouchEnv( q )
%t=TOUCHENV(q,n) When does SNR curve touch envelope?
%   q = vector of nearest neighbour transition rates, or length of chain

% if isscalar(q) && isint(q) && q>1
%     q=ones(1,q-1);
% end
% 
% [Wp,Wm,w]=MakeSMS(q);
[Wp,Wm,w]=DiffJump(q);

[qa,Ia]=SpectrumWpm(Wp,Wm,0.5,w);
Ia=Ia.*qa;

gamma=sqrt(128)/pi^2;

tm=(length(q)/3)^2;
options=optimset('Display','off');

tm=fminsearch(@MisMatch,tm,options);


    function s=SNR(t)
        s=Ia.'*exp(-qa*t);
    end

    function s=env(t)
        s=gamma/sqrt(2*exp(1)*t);
    end

    function ds=MisMatch(t)
        ds=env(t)-SNR(t);
    end


end

