function hess = Mats2ParamsHess( hesspp,hesspm,hessmp,hessmm )
%hess=MATS2PARAMSHESS( hesspp,hesspm,hessmp,hessmm )

n=size(hesspp,1);

hesspp=tens2mat(hesspp);
hesspm=tens2mat(hesspm);
hessmp=tens2mat(hessmp);
hessmm=tens2mat(hessmm);

hess=[hesspp hesspm; hessmp hessmm];

    function mat=tens2mat(tens)
        tens=permute(tens,[2 1 4 3]);
        mat=reshape(tens,n^2,[]);
        mat(1:n+1:end,:)=[];
        mat(:,1:n+1:end)=[];
    end

end