function [ mu, sigma ] = conditional_normal_sim(ICab,ind,Val)
mu=0;
sigma=0;
[nx,~]=size(ICab);
[n1,~]=size(Val);
if nx~=n1
    display('pas bonnes dimensions... pas content!! :-S');
    return;
end
var_krig=1/ICab(ind,ind);
mean_krig=0;

ICab2=[ICab(1:ind-1,ind);ICab(ind+1:end,ind)];
Val2=[Val(1:ind-1,1);Val(ind+1:end,1)];

Lab=-ICab2*var_krig;
mean_krig=Lab'*Val2;

mu=mean_krig;
sigma=sqrt(var_krig);


end

