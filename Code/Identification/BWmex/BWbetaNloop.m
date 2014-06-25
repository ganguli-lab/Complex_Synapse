function [ beta ] = BWbetaNloop( beta,updater )
%beta=BWBETANLOOP(eta,updater) normalised backward variables for Baum-Welch algorithm
%   beta     = normalised backward variables
%   updater  = updater matrices (M*outProj*eta)


siz=size(beta,1)*[1 1];

for i=size(beta,2):-1:2
    beta(:,i-1) = reshape(updater(:,:,i-1),siz) * beta(:,i);
end

end

