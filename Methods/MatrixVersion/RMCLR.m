function [L,S] = RMCLR(Y,Omega,lambda,rho,nu,MaxLoop,tol)
%RCLR 此处显示有关此函数的摘要
%   此处显示详细说明
%% Preparation
isize=size(Y);
L=Y;S=zeros(size(Y));
Z=L;
A=zeros(size(Y));B=zeros(size(Z));
cri=[];
%%  Main Loop
for loop=1:MaxLoop
    %Update L
    L=(Omega+1).^(-1).*(Omega.*(Y+A/rho-S)+Z+B/rho);
    %Update Z
    Aux=L-B/rho;
    [U,Sig,V]=svd(Aux,'econ');
    Z=U*diag(soft(diag(Sig),1/rho))*V';
    %Update S
    S=Omega.*soft(Omega.*(Y+A/rho-L),lambda/rho);
    
    %Multipliers
    A=A+rho*(Y-Omega.*(L+S));
    B=B+rho*(Z-L);
    rho=nu*rho;
    
    cri=[cri;norm(Y-Omega.*(L+S)),norm(Z-L)];
    if mod(loop, 20) == 0 | loop==1
        fprintf('Iteration No.%d，Convergence criterion=%f\n', loop, max(cri(end,:)));
    end
    if max(cri(end,:))<tol
        fprintf('======================Iteration No.%d，Algorithm converges======================\n', loop);
        break
    end
end
end

