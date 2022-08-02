function I2=fus_normalise(I1)

Q=size(I1,3);
if Q>1
    for i=1:Q
        I2(:,:,i)=I1(:,:,i)-min(min(I1(:,:,i)));
        I2(:,:,i)=I2(:,:,i)/max(max(I2(:,:,i)));    
    end
else
    I2=I1-min(min(I1));
    I2=I2/max(max(I2));
end

