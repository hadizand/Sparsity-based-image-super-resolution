clc;close all;clear
%% The mandrills eye
clear;clc;close all;
k=128;M=k; N=256; 
AI=imread('...\slice_10.PNG');

%IMG=load('mandrill'); AI=IMG.X;
 if size(AI,3)==3      Gs=rgb2gray(AI); else Gs=AI; end
IM_Original=(imresize(Gs,[N N]));
%IM_Original = imnoise(IM_Original,'gaussian',0,.00);

B=double(imresize(IM_Original,[k k]));
scale = N/(M);
ImResize_Nearest = imresize(B,scale,'nearest');
ImResize_bilinear = imresize(B,scale,'bilinear');
ImResize_bicubic = imresize(B,scale,'bicubic');

%%--------------------------------------CS
m=N/M;A=zeros(M,N);
for i=1:M %measurement matrix
A(i,1+(i-1)*m:(i)*m) = 1;
end

dict = wmpdictionary(N,'LstCpt',{'dct'});
%dict = wmpdictionary(N,'lstcpt',{{'db10',10}});
A1= A*dict;

for i=1:1 %initialization of recovery algorithm, smooth 0_norm (sl0)
    sigma_off = 0.001;
    A_pinv1 = pinv(A1);;mu_0 = 2;sigma_decrease_factor = 0.5;L = 3;
    true_s = sparseSigGen4plusNoise(9,floor(27/4),sigma_off);
    if sigma_off>0
        sigma_min = sigma_off*4;
    else
        sigma_min = 0.00001;
    end
end

for i=1:M%column interpolation
        xp=SL0(A1, m*double(B(i,:))', sigma_min, sigma_decrease_factor, mu_0, L, A_pinv1);
        %xp = l1eq_pd(A1'*double(B(i,:))', A1, [], m*double(B(i,:))', 1e-1);
        cl_recov(i,:)=dict*xp;       
end
for j=1:N
        xp=SL0(A1, m*double(cl_recov(:,j)), sigma_min, sigma_decrease_factor, mu_0, L, A_pinv1);
        %xp = l1eq_pd(A1'*double(zz1(:,j)), A1, [],m* double(zz1(:,j)), 1e-1);
        cs_rec(:,j)=dict*xp; 
        

end


%-------------------comparison----PSNR and MSE
peaksnr_CS = psnr(uint8(cs_rec),uint8(IM_Original));
peaksnr_Nearest = psnr(uint8(ImResize_Nearest),uint8(IM_Original));
peaksnr_bilinear = psnr(uint8(ImResize_bilinear),uint8(IM_Original));
peaksnr_bicubic = psnr(uint8(ImResize_bicubic),uint8(IM_Original));
display('[peaksnr_Nearest,peaksnr_bilinear,peaksnr_bicubic, peaksnr_CS]');
[peaksnr_Nearest,peaksnr_bilinear,peaksnr_bicubic, peaksnr_CS]
%---------MSE
immse_CS = immse(uint8(cs_rec),uint8(IM_Original));
immse_Nearest = immse(uint8(ImResize_Nearest),uint8(IM_Original));
immse_bilinear = immse(uint8(ImResize_bilinear),uint8(IM_Original));
immse_bicubic = immse(uint8(ImResize_bicubic),uint8(IM_Original));
display('[immse_Nearest,immse_bilinear,immse_bicubic, immse_CS]');
[immse_Nearest,immse_bilinear,immse_bicubic, immse_CS]

%-----------SSIM
ssim(uint8(cs_rec),uint8(IM_Original))
ssim_CS =ssim(uint8(cs_rec),uint8(IM_Original));
ssim_Nearest = ssim(uint8(ImResize_Nearest),uint8(IM_Original));
ssim_bilinear = ssim(uint8(ImResize_bilinear),uint8(IM_Original));
ssim_bicubic = ssim(uint8(ImResize_bicubic),uint8(IM_Original));

display('[ssim_Nearest,ssim_bilinear,ssim_bicubic, ssim_CS]');

[ssim_Nearest,ssim_bilinear,ssim_bicubic, ssim_CS]

%-------------------vidual comparison
figure;
h1=subplot(2,3,1);imshow(IM_Original,[]);title({'Original image ';['size: ',num2str(N),'*',num2str(N),'px']});
h2=subplot(2,3,2);imshow(B,[]);title({'Down-sampled ';['size: ', num2str(M),'*',num2str(M),'px']});linkaxes([h1,h2])
subplot(2,3,3);imshow(ImResize_Nearest,[]);title({['PSNR nearest: ', num2str(peaksnr_Nearest)];['MSE nearest: ', num2str(immse_Nearest)];['SSIM nearest: ', num2str(ssim_Nearest)]});
subplot(2,3,4);imshow(ImResize_bilinear,[]);title({['PSNR bilinear: ', num2str(peaksnr_bilinear)];['MSE bilinear: ', num2str(immse_bilinear)];['SSIM bilinear: ', num2str(ssim_bilinear)]});
subplot(2,3,5);imshow(ImResize_bicubic,[]);title({['PSNR bicubic: ', num2str(peaksnr_bicubic)];['MSE bicubic: ', num2str(immse_bicubic)];['SSIM bicubic: ', num2str(ssim_bicubic)]});
subplot(2,3,6);imshow(uint8(cs_rec),[]);title({['PSNR CS: ', num2str(peaksnr_CS)];['MSE CS: ', num2str(immse_CS)];['SSIM CS: ', num2str(ssim_CS)]})

sizes = [size(IM_Original);size(B); size(ImResize_Nearest);size(ImResize_bilinear);size(ImResize_bicubic);size(cs_rec)];

%--------------------to check how pixel to pixel CS and bicubic trace the
%initial image's pixcels
% figure;%whole pixels
% err_CS=double(IM_Original(:))-cs_rec(:);
% err_bicubic = double(IM_Original(:))-ImResize_bicubic(:);hold on
% plot(err_CS(:));plot(err_bicubic(:),'g');
figure;
x0=500;
y0=10;
width=550;
height=400;
set(gcf,'units','points','position',[x0,y0,width,height])

subplot(2,1,1);%Max of each column
plot(max(IM_Original));hold on
plot(max(cs_rec),'r');plot(max(ImResize_bicubic),'g');xlim([0 N]);
title('Max of each column');legend('original','bicubic','CS')

subplot(2,1,2);%Max of each row
plot(max(IM_Original'));hold on
plot(max(cs_rec'),'r');plot(max(ImResize_bicubic'),'g');xlim([0 N]);
title('Max of each row');legend('original','bicubic','CS')


%---------------------------------
ssim(uint8(cs_rec),uint8(IM_Original))

