clear all;close all;
addpath(genpath(pwd));
MRec_error=[];
MIR=[];
MOR=[];
MNu=[];
MP=[];
MSR=[];

for sam=10:-1:2
Rec_error=[];
IR=[];
OR=[];
Nu=[];
for i=1:20


N=500;
rank=13;
Le=randn(N,rank)/sqrt(N);Ri=randn(N,rank)/sqrt(N);
X=Le*Ri';
[U,~,V]=svds(X,rank);
alpha=1e-7;
Nu=[Nu,max(abs(U*V'),[],'all')*N*log(N)^alpha/sqrt(rank)/sqrt(N/rank)];

% % X=X./max(X(:));
Dim=size(X);

%%



ms = sam:sam:N-sam;
% % 
Omega = ones(N,N);
Omega(:,ms) = 0;
Omega(ms,:) = 0;
Omega=logical(Omega);




sample_ratio=sum(Omega(:))/N^2;
%% Noise supportive construction via golfing scheme
noise_ratio=0.05;
Cg=1;
k=floor(Cg*log(N));
Rho     = zeros(size(X));
D=-U*V';
conditionvalue=[];
for i=1:k
    conditionvalue=[conditionvalue,norm(P_T(U,V,~Omega.*(P_T(U,V,D))),'Inf')/norm(D,'Inf')];
    Rho_aux=zeros(size(X));
    eta_ratio=1-(noise_ratio)^(1/k);
    indd       = randperm(prod(size(X)));
    noised     = indd(1:floor(prod(size(X))*eta_ratio));
    Rho_aux(noised) = 1;
    D=P_T(U,V,D)-eta_ratio^(-1)*P_T(U,V,Omega.*Rho_aux.*P_T(U,V,D));
    Rho=Rho+Rho_aux;
end
IR=[IR,max(conditionvalue)]
%%
Rho    = ~(Rho > 0);
S=Rho.*randn(Dim);
Y=(X+S).*Omega;
indices=[];
%%
maxIter=1000;tol=1e-5;
[Eigen] = Pow_Iter(Omega,maxIter,tol,N,U,V);
OR=[OR,Eigen]
%% CO-RMCDS
lambda=1.5/sqrt((1-noise_ratio)*N)*log(N)^alpha;
rho=1e-2;
nu=1.1;
MaxLoop=5000;
tol=1e-5;
% % [L,S] = RCLR(Y,Omega,ksize,lambda,rho,nu,MaxLoop,tol);
[L,S] = RMCLR(Y,Omega,lambda,rho,nu,MaxLoop,tol);

Rec_error=[Rec_error,norm(L-X,'fro')/norm(X,'fro')]
end
MSR=[MSR,sample_ratio];
MRec_error=[MRec_error,mean(Rec_error)];
MIR=[MIR,mean(IR)];
MOR=[MOR,mean(OR)];
MNu=[MNu,mean(Nu)];
MP=[MP,mean(Rec_error<1e-4)];
end