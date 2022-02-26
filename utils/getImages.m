function [img_color  images] = getImages(img1 ,img2 )
if size(img1,3)==3
    img_color=cat(4,double(img1),double(img2));
    img1=rgb2gray(img1);
    img2=rgb2gray(img2);    
else
    img_color=cat(3,double(img1),double(img2));
end
    images=cat(3,double(img1),double(img2));
end

