function [u,v,Rtime] = optic_flow_compute(img1,img2,para1)
         [img_color images] = getImages(img1,img2 );
%% structure_texture_decomposition
         [images,images1] = structure_texture_decomposition_rof(images,1/8,100,para1.blend);    
    if para1.lambdaHOG~=0            
       [p,q,~]=size(img_color(:,:,:,1));
       para=get_para_flow(p,q);
       para.images=images;   
       [c1,c2,confidence] = getHOG(para );
       para.c1=c1;
       para.c2=c2;
       para.confidence=confidence;
    end
    para.lambdaH=para1.lambdaHOG;
%% Turn to LAB space
im1 = RGB2Lab(img_color(:,:,:,1));          
   for j = 1:size(im1, 3);
            im1(:,:,j) = scale_image(im1(:,:,j), 0, 255);
   end;     
        
      num_levels = para1.num_levels;
      smooth_sigma      = para1.smooth_sigma;  
      f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma) ;    
      pyramid_images    = compute_image_pyramid(images, f, num_levels , para1.space);
      pyramid_color_images    = compute_image_pyramid(im1, f, num_levels ,  para1.space);									
      u = zeros(size(img1,1),size(img1,2)); 
      v = zeros(size(img1,1),size(img1,2));
      uv(:,:,1)=u;
      uv(:,:,2)=v;
   if para1.nsigma
      para1.OCC=ones(size(pyramid_images{num_levels}, 1) ,size(pyramid_images{num_levels}, 2)); %% occlusion region
   end
   tic
 %% Flow field estiamtion using coarse-to-fine image pyramid warping technique
for i = num_levels: -1 : 1
    para1.i=i;
            nsz         = [size(pyramid_images{i}, 1) size(pyramid_images{i}, 2)];
            images      = pyramid_images{i};     
            img_color   = pyramid_color_images {i};  
            [u_descr, v_descr ,ro_descr,delta_inds] = getHOGUV(nsz(2),nsz(1),para );
    if i~=num_levels && para1.nsigma                 
       para1.OCC = ~imresize(OCC, nsz, 'nearest');
    end
            uv= resample_flow(uv, nsz);            
            u=uv(:,:,1);
            v=uv(:,:,2); 
           
       [ It Ikx2 Iky2 Ixz Iyz Ixx Ixy Iyy ] = getDerive( images, uv ,para1.h); 
       para1.Ixz=Ixz; para1.Iyz=Iyz; para1.Ixx=Ixx; para1.Ixy=Ixy; para1.Iyy=Iyy;
       [du, dv] = resolutionProcess_brox( It, Ikx2, Iky2, 1.8, u, v, 3, 500,u_descr,v_descr,ro_descr,delta_inds,para1.lambdaHOG,para1,images    ) ;
       uv(:,:,1)= uv(:,:,1)+du;
       uv(:,:,2)= uv(:,:,2)+dv;
  %% Improved median filter based on occlusion detection 
     if ~para1.nsigma
           [uv] = getMedianuv( uv ,  img_color ,images,para1);       
     else
           [uv OCC] = getMedianuv( uv ,  img_color ,images,para1);  
     end
         u=uv(:,:,1);
         v=uv(:,:,2);
       
end
 Rtime= toc/60;  
         