function [ Wdual ] = ReverseProcess( W,p )
%Wdual=REVERSEPROCESS(W,P) adjoint of W wrt measure defined by P
%   If W is a stochastic matrix and P is the equilibrium process, Wdual is
%   the time reversed process.
%   W = matrix/vector we're adjointing
%   P = measure in definition of adjoint (default = EqProb(W))

error(CheckSize(W,@ismat));

if isstochasticC(W) && ~exist('p','var')
    p=EqProb(W);
end

error(CheckValue(p,@isprob));
error(CheckValue(p,@(x) all(x>0),'non-degenerate'));
error(CheckSize(p,@(x) length(x)==length(W),'same length as W'));

Wdual=W';

if size(Wdual,1)>1
    Wdual = diag(1./p) * Wdual;
end
if size(Wdual,2)>1
    Wdual = Wdual * diag(p);
end


end

