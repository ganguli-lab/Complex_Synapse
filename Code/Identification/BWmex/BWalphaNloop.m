function [ alpha,eta,updater ] = BWalphaNloop( alpha,eta,updater )
%[alpha,eta]=BWALPHANLOOP(alpha,eta,updater) normalised forward variables for Baum-Welch algorithm
%   alpha    = normalised forward variables
%   eta      = normalisation factor
%   updater  = updater matrices (M*outProj*eta)

siz=size(alpha,2)*[1 1];

for i=2:length(eta)
    alpha(i,:) = alpha(i-1,:) * reshape(updater(:,:,i-1),siz);
    eta(i) = 1 / sum(alpha(i,:));
    alpha(i,:) = alpha(i,:) * eta(i);
    updater(:,:,i-1)=updater(:,:,i-1)*eta(i);
end
    
end

