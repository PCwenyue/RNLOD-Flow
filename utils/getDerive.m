function [ It Ikx2 Iky2 Ixz Iyz Ixx Ixy Iyy  ] = getDerive(images, uv ,h)
              [ It Ikx2 Iky2] = partial_deriv(images, uv);          
              I1x = imfilter(images(:,:,1), h,  'corr', 'symmetric', 'same');  
              I1y = imfilter(images(:,:,1),  h', 'corr', 'symmetric', 'same');
              Ixz=Ikx2-I1x;
              Iyz=Iky2-I1y;
              Ixx = 0.5*imfilter(Ikx2, h,  'corr', 'symmetric', 'same')+0.5*imfilter(I1x, h,  'corr', 'symmetric', 'same');
              Ixy = 0.5*imfilter(Ikx2, h', 'corr', 'symmetric', 'same')+0.5*imfilter(I1x, h',  'corr', 'symmetric', 'same');
              Iyy =0.5* imfilter(Iky2, h',  'corr', 'symmetric', 'same')+0.5*imfilter(I1y, h',  'corr', 'symmetric', 'same');
              
end

