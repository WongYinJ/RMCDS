function [Eigen] = Pow_Iter(Omega,maxIter,tol,N,U,V)
%POWER_ITER 此处显示有关此函数的摘要
%   此处显示详细说明
EigenM=rand(N);
for iter=1:maxIter
    Pre=EigenM;
    EigenM=P_T(U,V,~Omega.*(P_T(U,V,EigenM)));
    EigenM=EigenM/norm(EigenM,'fro');
    if norm(EigenM-Pre,'fro')<tol
        break
    end
end
iter
norm(EigenM-Pre,'fro')
Eigen=norm(P_T(U,V,~Omega.*(P_T(U,V,EigenM))),'fro');
end

