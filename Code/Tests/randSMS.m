
crq=zeros(1,100);
q=zeros(100,99);

w=[-ones(50,1);ones(50,1)];

for i=1:100
q(i,:)=rand(1,99);
[Wp,Wm]=MakeSMS(q(i,:));
[qa,ca]=Spectrum(Wp,Wm,0.5,w);
crq(i)=max(ca.*sqrt(qa));
end

clear w Wp Wm qa ca i;

[crq,ix]=sort(crq,'descend');
q=q(ix,:);