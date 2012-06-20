function [ h ] = TestCorder( len,fp,varargin )
%H=TESTCORDER(LEN,FP,...) test if area coeffs increase with random matrix
%   LEN = # states
%   FP = Fraction of potentiation transitions
%   H = handle
%   ... passed to bar

assert(isscalar(len));
assert(len>=0 && mod(len,2)==0);
assert(isscalar(fp));
assert(0<=fp && fp<=1);%fp in [0,1]

W=RandTrans(len);
w=ones(len,1);
w(1:(len/2))=-1;

[Wp,Wm,w]=TriangleDcmp(W,fp,w);

h=AreaCoeffBar(Wp,Wm,fp,w,varargin{:});

end

