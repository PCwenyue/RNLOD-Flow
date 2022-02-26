function  occusiondetection( uv,images)
addpath('utils');
  e1 = edge(uv(:,:,1), 'sobel');
  e2 = edge(uv(:,:,2), 'sobel');
  e  = e1|e2;
  mask = imdilate(e, ones([ 10 10]) );
  index = find(mask ==1);
  ht=size(uv,1);
            wt=size(uv,2);
            cc=zeros(ht,wt);
            Delaunay1=zeros(2*(ht-1)*(wt-1),15);
            count=1;
            [x,y]   = meshgrid(1:wt,1:ht);
            x2      = x + uv(:,:,1);        
            y2      = y + uv(:,:,2);  
            It  = partial_deriv(images, uv);
         for i=1:ht
            for j=1:wt
               if i~=ht && j~=wt
                  Delaunay1(count,1)=i; Delaunay1(count,2)=j;
              
                  Delaunay1(count,3)=i+ uv(i,j,2);Delaunay1(count,4)=i+uv(i,j+1,2);Delaunay1(count,5)=i+1+ uv(i+1,j,2);
                  Delaunay1(count,6)=j+uv(i,j,1);Delaunay1(count,7)=j+1+uv(i,j+1,1);Delaunay1(count,8)=j+uv(i+1,j,1);
                  Delaunay1(count,9)=abs(It(i,j));Delaunay1(count,10)=abs(It(i,j+1));Delaunay1(count,11)=abs(It(i+1,j));
                  Delaunay1(count,12)=i;  Delaunay1(count,13)=j+1;  Delaunay1(count,14)=i+1;  Delaunay1(count,15)=j;
                  count=count+1;
                   
                  Delaunay1(count,1)=i+1; Delaunay1(count,2)=j;
                  Delaunay1(count,3)=i+1+uv(i+1,j,2);Delaunay1(count,4)=i+1+uv(i+1,j+1,2);Delaunay1(count,5)=i+uv(i,j+1,2);
                  Delaunay1(count,6)=j+uv(i+1,j,1);Delaunay1(count,7)=j+1+uv(i+1,j+1,1);Delaunay1(count,8)=j+1+uv(i,j+1,1);
                  Delaunay1(count,9)=abs(It(i+1,j));Delaunay1(count,10)=abs(It(i+1,j+1));Delaunay1(count,11)=abs(It(i,j+1));
                  Delaunay1(count,12)=i+1;Delaunay1(count,13)=j+1;Delaunay1(count,14)=i;Delaunay1(count,15)=j+1;
                  count=count+1;
               end                
            end
        end   
   
        ll=size(Delaunay1,1);
        
        for f=1:ll
             if any((Delaunay1(f,1)+(Delaunay1(f,2)-1)*ht)==index)
       
             xr1=Delaunay1(f,1);yr1=Delaunay1(f,2);
    
            if xr1-25<1
                xl1=1;
            else
                xl1=xr1-25;
            end
            if xr1+25>ht
                xl2=ht;
            else
                xl2=xr1+25;
            end
            if yr1-25<1
                yl1=1;
            else
                yl1=yr1-25;
            end
            if yr1+25>wt
                yl2=wt;
            else
                yl2=yr1+25;
            end
            
            xc=[Delaunay1(f,3),Delaunay1(f,4),Delaunay1(f,5)];
            yc=[Delaunay1(f,6),Delaunay1(f,7),Delaunay1(f,8)];
            cf1=Delaunay1(f,9);
            cf2=Delaunay1(f,10);
            cf3=Delaunay1(f,11);
            [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
            [ii,jj]=find(IN-ON);
            if size(ii,1)==0
                continue;
            end
      
            for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                if cd<(cf1^2+cf2^2+cf3^2)/(cf1+cf2+cf3)
%                     
                    cc(Delaunay1(f,1),Delaunay1(f,2))=1;
                      cc(Delaunay1(f,12),Delaunay1(f,13))=1;
                       cc(Delaunay1(f,14),Delaunay1(f,15))=1;
                else
                    cc(xl1+i_-1,yl1+j_-1)=1;
                end
               
                
            end
                
            end
            
        end
         figure
 imshow(cc);

end

