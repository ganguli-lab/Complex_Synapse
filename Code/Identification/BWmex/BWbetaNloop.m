function [ beta ] = BWbetaNloop( beta,updater )
%beta=BWBETANLOOP(eta,updater) normalised backward variables for Baum-Welch algorithm
%   beta     = normalised backward variables
%   updater  = updater matrices (M*outProj*eta)

for i=size(beta,2):-1:2
    beta(:,i-1) = updater(:,:,i-1) * beta(:,i);
end

end

