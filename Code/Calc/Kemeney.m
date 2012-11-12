function [ eta ] = Kemeney( W )
%eta=KEMENEY(W) calculate Kemeney's constant
%   eta = sum_j T_ij p_j = Tr Z - tau 

error(CheckSize(W,@ismat));%matrix
error(CheckSize(W,@issquare));%square


qa=eig(-W);
qa=sort(qa);

eta = sum(1./qa(2:end));

end

