function [ binvec ] = BinVec( num,base,digits )
%bv=BINREP(num,base,digits) Vector of digits in Binary represantation
%   num = number to represent
%   base = base of representation (default=2)
%   digits = number of digits to use
%num =sum_i=1:digits bv(i)*base^(i-1)

existsAndDefault('base',2);
existsAndDefault('digits',max(floor(log(num)/log(base)),0)+1);

binvec=mod(floor(num./2.^(0:(digits-1))),2);

end

