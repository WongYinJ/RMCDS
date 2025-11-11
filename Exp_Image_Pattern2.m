clear all;close all;
addpath(genpath(pwd));

% 指定文件夹路径
folderPath = '.\RGB files';

% 获取文件夹中所有 .png 文件的列表
fileList = dir(fullfile(folderPath, '*.png'));

Indices1=[];Indices2=[];Indices3=[];Indices4=[];Indices5=[];Indices6=[];Indices11=[];Indices51=[];Indices61=[];
% 遍历文件列表
for k = 1:length(fileList)
    % 获取文件名
    fileName = fileList(k).name;
    
    % 构建完整文件路径
    filePath = fullfile(folderPath, fileName);
    
    % 读取图像
    img = imread(filePath);
    N =200;
    X=img(:,129:896,:);
    X=imresize(X,[N,N]);
    % 在这里对图像进行处理
    % 例如，显示图像
    X=double(rgb2gray(X));
    X=X/max(X(:));
%     imshow(X); 
    
    sam=5;
    ms = sam:sam:N-1;
    Omega = ones(N,N);
    Omega(:,ms) = 0;
    Omega(ms,:) = 0;
    Omega=logical(Omega);
    
    noise_ratio=0.01; %% rho
    X_ = imnoise(X,'salt & pepper',noise_ratio);

    Y=X_.*Omega;
    % % YY=X.*Omegaa;
%% CNNM
fprintf('======================CNNM======================\n');
ksize=[14,14];
L1 = inexact_alm_cnm_2D(Y, Omega,ksize);
Indices1=[Indices1;[psnr(L1,X),...
    rmse(L1,X),ssim(L1,X)]];
% Indices11=[Indices11;[psnr(double(uint8(L11)),X),rmse(double(uint8(L11)),X),ssim(double(uint8(L11)),X)]];
%% FGSR
fprintf('======================FGSR======================\n');
r=100;
options=[];
options.p=1/2^3;
options.d=ceil(r*1.5);
options.lambda=8e-1;
options.regul_B='L2';
[L2,~]=MC_FGSRp_PALM(Y,Omega,options);
Indices2=[Indices2;[psnr(L2,X),...
    rmse(L2,X),ssim(L2,X)]];
%% Lp-TNN
fprintf('======================Lp-TNN======================\n');
p = 1; r = 0; gamma = 1e-3; mu = 1e-1; rho = 1.01; maxIter = 50; tol = 1e-2;
L3 = Lp_TNN(Y, Omega, p, r, gamma, mu, rho, maxIter, tol);
Indices3=[Indices3;[psnr(L3,X),...
    rmse(L3,X),ssim(L3,X)]];
% imshow(L3,[]);
%% L0-BCD
fprintf('======================L0-BCD======================\n');
rak = 10;
maxiter = 8;
[L4, ~] = L0_BCD_F_image(X,Y,rak, maxiter);
Indices4=[Indices4;[psnr(L4,X),...
    rmse(L4,X),ssim(L4,X)]];
%% CO-RMCDS
fprintf('======================CO-RMCDS======================\n');
alpha=1e-7;
lambda=1/sqrt((1-noise_ratio)*N)*log(N)^alpha;
rho=1e-2;
nu=1.1;
MaxLoop=300;
tol=1e-5;
[L5,S5] = RMCLR(Y,Omega,lambda,rho,nu,MaxLoop,tol);
Indices5=[Indices5;[psnr(L5,X),...
    rmse(L5,X),ssim(L5,X)]];
% Indices51=[Indices51;[psnr(double(uint8(L51)),X),rmse(double(uint8(L51)),X),ssim(double(uint8(L51)),X)]];
%% GCO-RMCDS
fprintf('======================GCO-RMCDS======================\n');
ksize=[14,14];
% lambda=5e1/sqrt(N*log(N));
alpha=1e-7;
lambda=9/sqrt((1-noise_ratio)*N)*log(N)^alpha;
rho=1e-2;
nu=1.1;
mu=1.65e0;
MaxLoop=3000;
tol=1e-5;
[L6,S6] = RCLR_EquiADMM(Y,Omega,ksize,lambda,mu,rho,nu,MaxLoop,tol);
Indices6=[Indices6;psnr(L6,X),...
     rmse(L6,X),ssim(L6,X)];
end
MIndices1=mean(Indices1,1);
% MIndices11=mean(Indices11,1);
MIndices2=mean(Indices2,1);
MIndices3=mean(Indices3,1);
MIndices4=mean(Indices4,1);
MIndices5=mean(Indices5,1);
% MIndices51=mean(Indices51,1);
MIndices6=mean(Indices6,1);