function [ S ] = SNRcurveCaQa( t, ca,qa  )
%S=SNRCURVECaQa(T,ca,qa) SNR as function of time
%   T = time values
%   ca = area due to eigenmode
%   qa = decay rate of eigenmode

error(CheckSize(ca,@iscol));
error(CheckSize(qa,@(x) samesize(x,qa)));
error(CheckSize(t,@isrow));

% S = gmdmp(ca.*qa, 1, exp(-outer(qa,t,true)), 1);
S = (ca.*qa)'* exp(-qa*t);

end

