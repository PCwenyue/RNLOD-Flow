function [  u_descr, v_descr ,ro_descr,delta_inds] = getHOGUV(width,height,para )

if para.lambdaH==0
    u_descr=0; v_descr=0;ro_descr=0;delta_inds=0;
else
    c1= para.c1;c2=para.c2;confidence=para.confidence;
     centroids1_c=[c1(:,1)*width/para.q c1(:,2)*height/para.p];
    centroids2_c=[c2(:,1)*width/para.q c2(:,2)*height/para.p];
    c1r=round(centroids1_c);
    c2r=round(centroids2_c);
    keep=find(c1r(:,1)>0&c1r(:,2)>0&c1r(:,1)>0&c2r(:,2)>0 & ...
        c2r(:,2)<=height & c1r(:,2)<=height & c2r(:,1)<=width& c1r(:,1)<=width);
    centroids1_c=centroids1_c(keep,:);
    centroids2_c=centroids2_c(keep,:);
    confidence_c=confidence(keep);
    delta_inds=sub2ind([height width],round(centroids1_c(:,2)),...
        round(centroids1_c(:,1)));
    u_descr=zeros(width*height,1);
    u_descr(delta_inds)=centroids2_c(:,1)-centroids1_c(:,1);
    v_descr=zeros(width*height,1);
    v_descr(delta_inds)=centroids2_c(:,2)-centroids1_c(:,2);
    ro_descr=zeros(width*height,1);
    ro_descr(delta_inds)=confidence_c;
end
        
           
end

