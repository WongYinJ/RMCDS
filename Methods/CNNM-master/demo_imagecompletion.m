
%EXP3 此处显示有关此函数的摘要
%   此处显示详细说明
close all;
% % X0=load('E:\BaiduNetdiskWorkspace\BIT\Scripts\TensorCompletionBox\MTTD-main\test_data\test_colorimage\peppers.mat');
% % X0=load('testMSI_1.mat');

% % X=load('PaviaU.mat');
m =200;
% % X0=X.paviaU(1:m,1:m,10);

X0 = imread('lena.png');
X0 = imresize(X0,[m m]);

X0 = im2double(X0);
pval = max(X0(:));
X0=X0./pval;

ms = 5:5:195;

Omega = ones(m,m);
Omega(:,ms) = 0;
Omega(ms,:) = 0;
%%
Stripe=0.13*randn(m,m);
Ind=zeros(m,m);
ind=randperm(m);
Ind(:,ind(1:0.1*m))=1;
Stripe=Stripe.*Ind;

Gaussian=0.03*randn(m,m);

Sparse=0.03*randn(m,m);
sam_ratio=0.1;
Ome     = zeros(size(X0));
ind       = randperm(prod(size(X0)));
known     = ind(1:floor(prod(size(X0))*sam_ratio));
Ome(known) = 1;
Ome     = (Ome > 0);
Sparse=Sparse.*Ome;
%%


% % % % sample_ratio=0.1;
% % % % Omega     = zeros(size(X0));
% % % % ind       = randperm(prod(size(X0)));
% % % % known     = ind(1:floor(prod(size(X0))*sample_ratio));
% % % % Omega(known) = 1;
% % % % Omega     = (Omega > 0);

% % % Omega=logical(CFA_Mask([100,100],[2,2]));


inds = Omega>0.1;
disp(['missing rate = ' num2str(sum(~inds(:))/(m*m*9))]);
Y = X0.*Omega;
subplot(1,3,1);
imshow(Y,[]);
title('input');

%run DFT_L1
% % % tic;
% % % X_hat = inexact_alm_dft_nD(Y,Omega);
% % % toc;
% % % subplot(1,3,2);
% % % imshow(X_hat(:,:,1));
% % % title('DFT_L1');
% % % % % acc = psnr(X_hat(~inds),X0(~inds),pval);
% % % % % disp(['PSNR (DFT_L1) = ' num2str(acc)]);


% run_CNNM
ksize = [13 13,4];
tic;
X_hat = inexact_alm_cnm_3D(Y, Omega,ksize);
toc;
subplot(1,3,3);
imshow(X_hat(:,:,1));
title('CNNM');
% % acc = psnr(X_hat(~inds),X0(~inds),pval);
% % disp(['PSNR (CNNM) = ' num2str(acc)]);
rmse=sqrt(mean((X_hat(:)-X0(:)).^2))

