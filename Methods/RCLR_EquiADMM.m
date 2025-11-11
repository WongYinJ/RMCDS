function [L,S] = RCLR_EquiADMM(Y,Omega,ksize,lambda,mu,rho,nu,MaxLoop,tol)
%RCLR 此处显示有关此函数的摘要
%   此处显示详细说明
%% Preparation
isize=size(Y);
L=Y;S=zeros(size(Y));
Z=cconv2mtx(L,ksize);
B=zeros(size(Z));
cri=[];ps=[];
%%  Main Loop
for loop=1:MaxLoop
    %Update L
    L=(mu*Omega+rho).^(-1).*(mu*Omega.*(Y-S)+adj_2D(rho*Z+B, isize, ksize)/prod(ksize));
    %Update Z
    Aux=cconv2mtx(L,ksize)-B/rho;
    [U,Sig,V]=svd(Aux,'econ');
    Z=U*diag(soft(diag(Sig),1/rho))*V';
    %Update S
    S=Omega.*soft(Omega.*(Y-L),lambda/mu/prod(ksize));
    
    %Multipliers
    B=B+rho*(Z-cconv2mtx(L,ksize));
    rho=nu*rho;
    
    cri=[cri;norm(Z-cconv2mtx(L,ksize))];

    
    if mod(loop, 20) == 0 | loop==1
        fprintf('Iteration No.%d，Convergence criterion=%f\n', loop, cri(end));
    end
    if cri(end)<tol
        fprintf('======================Iteration No.%d，Algorithm converges======================\n', loop);
        break
    end
end
end

