%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Demo file to demonstrate the ICA fusion scheme 
%%%%%%%% This demo assumes that ICA bases have been trained,
%%%%%%%% i.e. the icabases.mat file exists. If not, run
%%%%%%%% the ICA_bases_train.m file first.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Nikolaos Mitianoudis
%%%%%%% Democritus University of Thrace, 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
%%%%%%%%%% Load images to be fused %%%%%
I1=double(imread('./fusion_images/Disk1.gif'))./255; 
I2=double(imread('./fusion_images/Disk2.gif'))./255 ;
x(:,:,1)=fus_normalise(I1); %% Stretch image contrast
x(:,:,2)=fus_normalise(I2);
figure
subplot(2,1,1);imagesc(x(:,:,1));colormap(gray);axis equal;axis off;
title(' Input 1')
subplot(2,1,2);imagesc(x(:,:,2));colormap(gray);axis equal;axis off;
title(' Input 2');


load icabases;

%%%%%%%%%% Simple fusion examples with different rules.
%%%%%%%%% Fuse the two images using the "maxabs" rule, without optimal
%%%%%%%%% contrast and no denoising. 
u1=ica_fusion_prog_mn(x,frame,Picao,'maxabs',0,0,0);

figure; 
 imagesc(u1);colormap(gray);axis equal;axis off;
title('maxabs-no optimal contrast- no denoising');
pause;
%%%%%%%%% Fuse the two images using the "regional" rule and the "entropy"
%%%%%%%%% texture features without optimal contrast and no denoising. 
u2=ica_fusion_prog_mn_june13(x,frame,Picao,'entr',0,0,0);

figure; 
imagesc(u2);colormap(gray);axis equal;axis off;
title('regional-entropy feat.-no optimal contrast- no denoising');
pause;

%%%%%%%%%% Simple fusion examples with denoising.
%%%% Add AGWN to images
[M N Q]=size(x);
x=x+0.05*randn(M,N,Q);
figure
subplot(2,1,1);imagesc(x(:,:,1));colormap(gray);axis equal;axis off;
title('Noisy Input 1')
subplot(2,1,2);imagesc(x(:,:,2));colormap(gray);axis equal;axis off;
title('Noisy Input 2');
%%%%%%%%% Fuse the two images using the "maxabs" rule, without optimal
%%%%%%%%% contrast and denoising with hard threshold 0.2 . 
u3=ica_fusion_prog_mn(x,frame,Picao,'maxabs',0,1,0.2);
figure; 
 imagesc(u3);colormap(gray);axis equal;axis off;
title('maxabs-no optimal contrast- with denoising');
pause;

%%%%%%%%% Multi-modal Example with optimal contrast
clear x;

I1=double(imread('./fusion_images/1817i.bmp'))./255; 
I2=double(imread('./fusion_images/1817v.bmp'))./255 ;
x(:,:,1)=fus_normalise(I1);
x(:,:,2)=fus_normalise(I2);


figure
subplot(2,1,1);imagesc(x(:,:,1));colormap(gray);axis equal;axis off;
title('Multimodal Input 1')
subplot(2,1,2);imagesc(x(:,:,2));colormap(gray);axis equal;axis off;
title('Multimodal Input 2');

%%%%%%%%% Fuse the two images using the "maxabs" rule, without optimal
%%%%%%%%% contrast and no denoising. 
u4=ica_fusion_prog_mn(x,frame,Picao,'maxabs',0,0,0);
figure; 
imagesc(u4);colormap(gray);axis equal;axis off;
title('maxabs - no optimal contrast - no denoising');
pause;
%%%%%%%%% Fuse the two images using the "maxabs" rule, with optimal
%%%%%%%%% contrast and no denoising. 
u4=ica_fusion_prog_mn(x,frame,Picao,'maxabs',1,0,0);
figure; 
imagesc(u4);colormap(gray);axis equal;axis off;
title('maxabs - with optimal contrast - no denoising');
