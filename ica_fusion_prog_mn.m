function  u1=ICA_fusion_prog_mn(x,frame,Picao,method,opt_contr,denoise,thres)

%%%%%%%% This the image fusion routine using ICA bases
%%%%%%% INPUTS: 
%%%%%%%         x   are the input images to be used arranged in a 3D array.
%%%%%%%             x(:,:,i) represents the i-th grayscale image    
%%%%%%%     frame   patch frame size (frame x frame)
%%%%%%%     Picao   Forward ICA transformation matrix, as estimated by
%%%%%%%             the ICAbases function.
%%%%%%%    method   defines the fusion method:
%%%%%%%                 'maxabs'            Max_abs fusion method
%%%%%%%                 'mean'              Mean fusion method
%%%%%%%                 'weighted'          Weighted Combination fusion method
%%%%%%%                 'regional'          Regional fusion method
%%%%%%%                 'adaptive'          Adaptive fusion method
%%%%%%%   opt_contr  If set to 1, it performes automated contrast
%%%%%%%  
%%%%%%%    denoise   if set to 1, it turns ON the denoising capability by
%%%%%%%              hard thresholding
%%%%%%%     thres    Has to be determined manually, from a constant area of
%%%%%%%              the image. It represents the variance of noise in the
%%%%%%%              ICA domain.
%%%%%%% OUTPUTS:
%%%%%%%    u1        the final fused image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%        Nikolaos Mitianoudis and Tania Stathaki 
%%%%%%%         Imperial College London  2008
%%%%%%%
%%%%%%% Please refer to the following papers, when using the software:
%%%%%%% [1] Mitianoudis N., Stathaki T., " Pixel-based and Region-based 
%%%%%%%     Image Fusion schemes using ICA bases", Information Fusion 8 (2), 
%%%%%%%     pp. 131-142, April 2007.
%%%%%%% [2] Mitianoudis N., Stathaki T., "Optimal Contrast Correction for 
%%%%%%%     ICA-based Fusion of Multimodal Images" , IEEE Sensors Journal, 
%%%%%%%     Vol. 8, No. 12, pp. 2016 - 2026, Dec. 2008. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N,M,Q]=size(x);
x_aug=zeros(N+frame*2-1,M+frame*2-1,Q);
x_aug(frame:frame+N-1,frame:frame+M-1,:)=x;
x=[];
x=x_aug;x_aug=[];[N,M,Q]=size(x);

disp('Segmenting images into patches');
LMLV=[];LV=[];
P=zeros(size(Picao,1),(N-frame+1)*(M-frame+1),Q);
for i=1:Q
  disp([ 'Image ' int2str(i)]);
  X1=im2col(x(:,:,i),[frame frame],'sliding');
  LMLV=[LMLV;mean(X1)];
  X1=X1-ones(frame*frame,1)*LMLV(i,:);
  LV=[LV;var(X1)];
  P(:,:,i)=Picao*X1;
end

%clear X1;




if denoise==1 
       disp('Denoising ON')
       P(find(abs(P)<thres))=0;
   else
       disp('Denoising OFF')
end


switch lower(method)
    case {'maxabs'}
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
        clear P; P=P_i; clear P_i; clear tmp;
        
    case 'mean'
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
        clear P; P=P_fin; clear P_fin;
    case 'weighted'
        disp('Weighted Combination fusion method');
        for i=1:size(P,3)
            E(i,:)=sum(abs(P(:,:,i)));
        end
        E_all=sum(E,1);
        E_all(find(E_all<0.0001))=0.0001;
        E=E./(ones(size(E,1),1)*E_all);
        clear E_all;
        P_fin=zeros(size(P,1),size(P,2));
        for i=1:size(P,3)
           P_fin=P_fin+P(:,:,i).*(ones(size(P,1),1)*E(i,:));         
        end
        if opt_contr==1
        Lsigf=sum((pinv(Picao)*P_fin).^2);
        Lsigxf=zeros(Q,length(Lsigf));
        for i=1:Q
            tmp=pinv(Picao)*squeeze(P(:,:,i));
            Lsigxf(i,:)=sum(abs(tmp.*(pinv(Picao)*P_fin)))./max(sum(tmp.^2)+Lsigf,0.0001);
        end
        end
        clear P; P=P_fin; clear P_fin;
        
    case 'regional'
        disp('Regional fusion method');
        MSK=zeros(size(P,3),size(P,2));
        
        for i=1:size(P,3)
            E(i,:)=sum(abs(P(:,:,i)));
            me(i)=2*mean(E(i,:));
            MSK(i,find(E(i,:)>me(i)))=1;
        end
        MSK=sum(MSK);MSK(find(MSK>0))=1;
        t1=find(MSK==1);
        t2=find(MSK==0);
        E_all=sum(E,1);
        E_all(find(E_all<0.0001))=0.0001;
        E=E./(ones(size(E,1),1)*E_all);
        clear E_all;
         P_fin=zeros(size(P,1),size(P,2));
        % P_fin=mean(P,3);
         for i=1:size(P,3)
           P_fin(:,t1)=P_fin(:,t1)+P(:,t1,i).*(ones(size(P,1),1)*E(i,t1));         
        end
         P_fin(:,t2)=mean(P(:,t2,:),3) ;
         
        if opt_contr==1
        Lsigf=sum((pinv(Picao)*P_fin).^2);
        Lsigxf=zeros(Q,length(Lsigf));
        for i=1:Q
            tmp=pinv(Picao)*squeeze(P(:,:,i));
            Lsigxf(i,:)=sum(abs(tmp.*(pinv(Picao)*P_fin)))./max(sum(tmp.^2)+Lsigf,0.0001);
        end
        end
         clear P; P=P_fin; clear P_fin;
       
        
    case 'adaptive'
        disp('Adaptive fusion method')
        P1=zeros(size(Picao,1),(N-frame+1)*(M-frame+1));
        for i=1:size(P,2)
            x=squeeze(P(:,i,:));x=x';
            iter=0;
            w=ones(size(x,1),1)/size(x,1);
            ep=ones(size(x,1),1);
            lambda=0.1;
            while iter <100
                iter=iter+1;
                y=w'*x;
                dj=sign(y)*x';        
                dj=(dj')/size(P,1);
                w=w+0.05*dj;
                w=abs(w)./sum(abs(w));
            end
            P1(:,i)=y';
        end
         if opt_contr==1
        Lsigf=sum((pinv(Picao)*P1).^2);
        Lsigxf=zeros(Q,length(Lsigf));
        for i=1:Q
            tmp=pinv(Picao)*squeeze(P(:,:,i));
            Lsigxf(i,:)=sum(abs(tmp.*(pinv(Picao)*P1)))./max(sum(tmp.^2)+Lsigf,0.0001);
        end
        end
        
        P=P1;clear P1;
    otherwise
        disp('Unknown method.')
end


U1=pinv(Picao)*P;




clear P;
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