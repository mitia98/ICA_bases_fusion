function  u1=ICA_fusion_prog_mn_june13(x,frame,Picao,region,opt_contr,denoise,thres)

%%%%%%%% This the image fusion routine using ICA bases
%%%%%%% INPUTS: 
%%%%%%%         x   are the input images to be used arranged in a 3D array.
%%%%%%%             x(:,:,i) represents the i-th image    
%%%%%%%     frame   patch frame size (frame x frame)
%%%%%%%     Picao   Forward ICA transformation matrix, as estimated by
%%%%%%%             the ICAbases function.
%%%%%%%     region   Selects the local texture information to be extracted,
%%%%%%%                 ('entropy' or 'std' or 'fft')
%%%%%%%   opt_contr  If set to 1, it performes automated contrast
%%%%%%%  
%%%%%%%    denoise   if set to 1, it turns ON the denoising capability by
%%%%%%%              hard thresholding
%%%%%%%     thres    Has to be determined manually, from a constant area of
%%%%%%%              the image. It represents the variance of noise in the
%%%%%%%              ICA domain.
%%%%%%% OUTPUTS:
%%%%%%%    u1        the final fused image

pp1=0.4; %%% Power of L1-norm for the segmentation activity detector
[N,M,Q]=size(x);
x_aug=zeros(N+frame*2-1,M+frame*2-1,Q);
x_aug(frame:frame+N-1,frame:frame+M-1,:)=x;
x=[];
x=x_aug;x_aug=[];[N,M,Q]=size(x);
iptsetpref('ImshowBorder','tight')
Wentr=zeros(Q,(N-frame+1)*(M-frame+1));
disp('Segmenting images into patches');
LMLV=[];LV=[];
P=zeros(size(Picao,1),(N-frame+1)*(M-frame+1),Q);
for i=1:Q
  disp([ 'Image ' int2str(i)]);
  
  switch lower(region)
          case {'entropy','entr'}
            tmp=entropyfilt(x(:,:,i));
          case 'std'
            tmp=stdfilt(x(:,:,i));
          case {'fft','fourier'}
             tmp=fourier_filt(x(:,:,i));
          otherwise
            disp('Unknown regional method')
   end
  

  tmp=tmp(1:(N-frame+1),1:(M-frame+1));
  Wentr(i,:)=reshape(tmp,1,(N-frame+1)*(M-frame+1));
  X1=im2col(x(:,:,i),[frame frame],'sliding');
  LMLV=[LMLV;mean(X1)];
  X1=X1-ones(frame*frame,1)*LMLV(i,:);
  LV=[LV;var(X1)];
  P(:,:,i)=Picao*X1;
end

Wentr_m=max(sum(Wentr),0.00001);
Wentr=Wentr./(ones(Q,1)*Wentr_m);
clear Wentr_m; clear  tmp;

%clear X1;

L1=zeros(Q,size(P,2));
for i=1:Q
    L1(i,:)=sum(abs(P(:,:,i)).^pp1);
    tmp2=L1(i,:);
 [tmp, ctrs] = kmeans(L1(i,:),3);
% [ctrs,U,obj_fcn] = fcm(L1(i,:)',3);

[~, idx]=sort(ctrs,1,'ascend');
 cl=[1 2 3]/3;
 
%  maxU=max(U);
 for j=1:3
     
      t=find(tmp==idx(j));
      L1(i,t)=j;
%       t=find(U(j,:) == maxU);
%      L1(i,t)=idx(j);
 end   
% figure;[p,q]=hist(tmp2,100); plot(q,p);
% figure;imagesc(reshape(L1(i,:),N-frame+1,M-frame+1));colormap(gray);axis equal;axis off;
end
L1m=max(L1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;[p,q]=hist(L1m,100); plot(q,p);
% [tmp, ctrs] = kmeans(L1m,3);
% 
% [~, idx]=sort(ctrs,1,'ascend');
% 
% cl=[1 2 3]/3;
% 
%  for j=1:3
%     t=find(tmp==idx(j));
%     L1m2(t)=j;
% end   

%figure;imshow(reshape(L1m,N-frame+1,M-frame+1)/3);%colormap(gray);axis equal;axis off;

t1= L1m==1;
t2=(L1m==2);
t3=(L1m==3);


ct=[];tmp=[];L1=[];

if denoise==1 
       disp('Denoising ON')
       P(find(abs(P)<thres))=0;
   else
       disp('Denoising OFF')
end


P_fus=zeros(size(P,1),size(P,2));


        disp('Max_abs fusion method')
        [P_max I]=max(abs(P),[],3); 
        clear P_max;
        P_i=zeros(size(I));
        for i=1:size(I,1)
            for j=1:size(I,2)
                P_i(i,j)=P(i,j,I(i,j));
            end
        end
        if opt_contr==1
        Lsigf=sum((pinv(Picao)*P_i).^2);
        Lsigxf=zeros(Q,length(Lsigf));
        for i=1:Q
            tmp=pinv(Picao)*squeeze(P(:,:,i));
            Lsigxf(i,:)=(sum(tmp.*(pinv(Picao)*P_i)))./max(sum(tmp.^2)+Lsigf,0.0001);
            
        end
        end
         P_fus(:,t3)=P_i(:,t3); clear P_i; clear tmp; 
        
    
        disp('Mean fusion method')
        P_fin=mean(P,3);
         if opt_contr==1
        Lsigf=sum((pinv(Picao)*P_fin).^2);
        Lsigxf=zeros(Q,length(Lsigf));
        for i=1:Q
            tmp=pinv(Picao)*squeeze(P(:,:,i));
            Lsigxf(i,:)=sum(abs(tmp.*(pinv(Picao)*P_fin)))./max(sum(tmp.^2)+Lsigf,0.0001);
        end
        end
         P_fus(:,t1)=P_fin(:,t1); clear P_fin;
    
        disp(['Regional method' region]);
%         for i=1:size(P,3)
%             E(i,:)=sum((P(:,:,i)).^2);
%         end
%         E_all=sum(E,1);
%         E_all(find(E_all<0.0001))=0.0001;
%         E=E./(ones(size(E,1),1)*E_all);
        clear E_all;
        P_fin=zeros(size(P,1),size(P,2));
        for i=1:size(P,3)
           P_fin=P_fin+P(:,:,i).*(ones(size(P,1),1)*Wentr(i,:));         
        end
        if opt_contr==1
        Lsigf=sum((pinv(Picao)*P_fin).^2);
        Lsigxf=zeros(Q,length(Lsigf));
        for i=1:Q
            tmp=pinv(Picao)*squeeze(P(:,:,i));
            Lsigxf(i,:)=sum(abs(tmp.*(pinv(Picao)*P_fin)))./max(sum(tmp.^2)+Lsigf,0.0001);
        end
        end
        P_fus(:,t2)=P_fin(:,t2); clear P_fin; clear P;
        
   


U1=pinv(Picao)*P_fus;

clear P_fus;
w=ones(1,Q)/Q;

if opt_contr==1
LV_m=sum(LV); LV_m(find(LV_m<0.00001))=0.00001;
LV=LV./(ones(Q,1)*LV_m);
%LV_m=max(LV); %LV_m=LV_m/sum(LV_m);
LV=LV.*Lsigxf;

DW=[];
  for i=1:100
    i
    mf=w*LMLV;
    tmp=0;
    for i=1:Q
      tmp=tmp+2*LV(i,:).*LMLV(i,:).*(LMLV(i,:).^2-mf.^2)./(max(LMLV(i,:).^2+mf.^2,0.00001).^2);
%      tmp=tmp+LMLV(i,:).*(LMLV(i,:).^2-mf.^2)./max(LMLV(i,:).^2+mf.^2,0.00001);
    end
    %tmp=LV_m.*tmp;
    tmp=LMLV.*(ones(Q,1)*tmp);
    dw=mean(tmp'); DW=[DW sum(dw.^2)]; 

    
    w=w+2.5*(dw);
    w=abs(w)./sum(abs(w))
    
  end
%plot(DW);pause(0.5);
else
    mf=w*LMLV;
end

w
LMLV=mf;
%%%%%%%%%%%%%%%%%  Reconstruct the fused image  
u1=zeros(N,M); 
for j=1:M-frame+1
    for i=1:N-frame+1
        p=U1(:,(j-1)*(N-frame+1)+i);         
        u1(i:i+frame-1,j:j+frame-1)=u1(i:i+frame-1,j:j+frame-1)+reshape(p,frame,frame)+LMLV((j-1)*(N-frame+1)+i);%*LMLV1(2,(i-1)*(M-frame)+j)
    end
end
u1=u1/(frame^2);
 u1(N-frame+1:N,:)=[]; u1(:,M-frame+1:M)=[];
 u1(1:frame-1,:)=[];u1(:,1:frame-1)=[];








% % 