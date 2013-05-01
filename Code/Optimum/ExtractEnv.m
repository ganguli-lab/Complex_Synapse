function [env,whichn,whicht,new_res,snr] = ExtractEnv( res_struct )
%EXTRACTENV(res_struct) extract envelope from grid search results struct
%   Detailed explanation goes here

snr=[res_struct.SNR];
snr=reshape(snr,[length(res_struct(1).SNR) size(res_struct) ]);
snr=real(snr);

env=snr;
[env,whicht]=max(env,[],3);
[env,whichn]=max(env,[],2);
whicht=diag(whicht(:,whichn));


new_res=res_struct(sub2ind(size(res_struct),whichn,whicht));


end

