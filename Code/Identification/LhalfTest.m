function [ val ] = LhalfTest( gamma,Lrow )
%LHALFTEST Summary of this function goes here
%   Detailed explanation goes here

        %correct value of Lagrange multiplier given by root of this function
        %   Lrow=Mrow/beta^2
        %   gamma^2 = 1/lambda
        onesrow=ones(size(Lrow));
        onescol=ones(size(gamma));

        %   Lnew=Mnew/beta^2
        %   Lrow=Mrow/beta^2
        Lnew = ( 2 * gamma.^2 * Lrow ...
            + gamma.^4 * onesrow ...
            - ( gamma.^3 * onesrow ) ...
            .* sqrt( gamma.^2 * onesrow + 4 * onescol * Lrow ) ) / 2 ;

        val = sum( Lnew, 2 ) - 1 ;


end

