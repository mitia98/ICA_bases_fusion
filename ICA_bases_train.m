%%%%%%%     Train ICA bases from natural images
%%%%%%%     Put all training images in the folder "training_images"
%%%%%%%     and run this batch file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Nikolaos Mitianoudis
%%%%%%  Democritus University of Thrace, 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Variables that can be adjusted.
frame=8;            %%%%% Size of extracted patches frame x frame
Nco=60;             %%%%% Number of ICA bases extracted. Has to be less than frame^2.
NPatches=1000;      %%%%% Number of extracted frames per image



%%%%%%%%%%% Main programme 
files = dir(fullfile('.\training_images\*.jpg'));
NoI = length(files);  %%%  Number of total input images used for training

P=zeros(frame^2,NoI*NPatches);


for k=1:NoI

 I1=double(imread(['.\training_images\'  files(k).name]))./255 ;
 if size(I1,3)~=1
    I1=rgb2gray(I1);
 end
[n1,n2]=size(I1);
I1=double(I1);

%%%% Randomly select patches from each image, normalise to zero mean and
%%%% perform lexicographic ordering, i.e. arrange in a vector.

for i=1:NPatches
    x=[floor(random('unif',1,n1-frame,1,1));floor(random('unif',1,n2-frame,1,1))];
    ppatch=I1(x(1):x(1)+frame-1,x(2):x(2)+frame-1);
   
    ppatch=reshape(ppatch,1,frame*frame);
     ss=sqrt(var(ppatch));
    if ss>0.0001
            ppatch=(ppatch-mean(ppatch));%./ss;
            P(:,NPatches*(k-1)+i)= ppatch';
      end
end

end

%%%%% Finally P is matrix containing all the extracted patches from the
%%%%% input images. Each column of P represents a different patch.
Picao=icabases(P,Nco,frame);
 savefile = 'icabases.mat';
 save(savefile, 'frame', 'Picao');
 disp('ICA bases file saved!!!')
