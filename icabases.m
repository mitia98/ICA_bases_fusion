
function Picao=ICAbases(P,Nco,frame)
%%%%%% This function takes the inputs from the extract_patches function 
%%%%%  and more especially the matrix P and trains ICA bases
%%%%%  INPUTS:   P     the matrix containing the extracted patches (patch/column)
%%%%%          frame   patch size frame x frame, as used in the extract_patches function 
%%%%%           Nco     Number of desired ICA bases (maximum choice is frame^2-1)
%%%%% OUTPUTS: Picao   the matrix size Nco x frame^2  that contains the forward ICA transformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Nikolaos Mitianoudis and Tania Stathaki
%%%%% Imperial College London 2008.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('PCA step.........');
%% PCA algorithm, first dimensionality reduction
[NN dsa]=size(P);
 MP=mean(P')';
 P=P-MP*ones(1,dsa);
 C=P*P'/dsa;
 [V,D]=svd(C);
 D1=diag(D);
 clear C;
 


Vm=V*real(diag(1./max(sqrt(D1),0.000001)));
V=[];
V=Vm(:,1:Nco);
Pred=(V')*P;
  
disp('ICA step.........');
  
%%%% ICA algorithm independence introduced

W=[];
b=zeros(Nco,1);
h = waitbar(0,'Estimating ICA bases...');
for i=1:Nco
     waitbar(i/Nco);
   w=2*rand(Nco,1)-1; 
   w=w/norm(w);    	      	 % Initial w that has norm(w)=1
   flag=0;
   iter=0;
   while flag==0 & iter <400
      u=w'*Pred;
      w1=Pred*tanh(u)'/dsa-mean(1-tanh(u).^2)*w;        

      
      w1=w1 -b*(b')*w1;
      w1 = w1/norm(w1);
      if abs((w1')*w) >=0.99999
         flag=1;
      end
      w=w1;
      iter=iter+1;
   end
   if i==1
      b=w;
   else  
      b=[b w];
   end
   
   W=[W ;w'];   
end
 W=W*real((W'*W)^(-0.5));
 Picao=W*(V');
 Picao=Picao*real((Picao'*Picao)^(-0.5));
 close(h)

 
 
 
 
     
  
