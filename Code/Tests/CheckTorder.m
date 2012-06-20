function [ tf,failures ] = CheckTorder( T )
%TF=CHECKTORDER(T) checks if T_ik > T_jk when eta^+(i>j>k)
%   assumes states are in order of decreasing eta^+

assert(ismat(T));
assert(issquare(T));
n=length(T);

tf=true(n,n,n);

for i=1:n
    for j=1:n
        for k=1:n
            if i<j
                if k>j
                    tf(i,j,k)=T(i,k)>T(j,k);
                elseif k<i
                    tf(i,j,k)=T(i,k)<T(j,k);
                end%if k
            elseif i>j
                if k<j
                    tf(i,j,k)=T(i,k)>T(j,k);
                elseif k>i
                    tf(i,j,k)=T(i,k)<T(j,k);
                end%if k
            end%if i
        end%for k
    end%for j
end%for i

[i,j,k]=ind2sub(size(tf),find(~tf));
failures=[i,j,k];

end

