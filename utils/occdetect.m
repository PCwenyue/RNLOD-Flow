function [cc] = occdetect(im1,im2,uv ,method,oc )
addpath('utils');

im1=rgb2gray(im1);
im2=rgb2gray(im2);
image(:,:,1)=double(im1);
image(:,:,2)=double(im2);

e1 = edge(uv(:,:,1), 'sobel');
e2 = edge(uv(:,:,2), 'sobel');
e  = e1|e2;
mask = imdilate(e, ones([ 10 10]) );
index = find(mask ==1);

if nargin == 5
  countf=find(oc(index)==255);
  countf=length(countf);
end
  ht=size(uv,1);
             wt=size(uv,2);
            cc=zeros(ht,wt);
            Delaunay1=zeros(2*(ht-1)*(wt-1),15);
             count=1;
            [x,y]   = meshgrid(1:wt,1:ht);
             x2      = x + uv(:,:,1);        
             y2      = y + uv(:,:,2);  
         
           
            It  = partial_deriv(image, uv);
            It_left=interp2(x,y,It,x+0.5,y,method ); 
            It_bottom=interp2(x,y,It,x,y+0.5,method ); 
            It_middle=interp2(x,y,It,x+0.5,y+0.5,method ); 
           u_left= interp2(x,y,uv(:,:,1),x+0.5,y,method); 
           v_left= interp2(x,y,uv(:,:,2),x+0.5,y,method); 
           u_bottom= interp2(x,y,uv(:,:,1),x,y+0.5,method ); 
           v_bottom= interp2(x,y,uv(:,:,2),x,y+0.5,method );
           u_middle= interp2(x,y,uv(:,:,1),x+0.5,y+0.5,method ); 
           v_middle= interp2(x,y,uv(:,:,2),x+0.5,y+0.5,method);
%            'nearest' - nearest neighbor interpolation
%     'linear'  - bilinear interpolation
%     'spline'  - spline interpolation
%     'cubic'   - bicubic interpolation
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
                    if rem(f,2)==1
                %1
                xc=[Delaunay1(f,3),Delaunay1(f,1)+0.5+v_bottom(Delaunay1(f,1),...
                Delaunay1(f,2)),Delaunay1(f,1)+v_left(Delaunay1(f,1),Delaunay1(f,2))];
                yc=[Delaunay1(f,6),Delaunay1(f,2)+u_bottom(Delaunay1(f,1),...
                Delaunay1(f,2)),Delaunay1(f,2)+0.5+u_left(Delaunay1(f,1),Delaunay1(f,2))];
                [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                [ii,jj]=find(IN-ON);
                 if size(ii,1)==0
            
                 
                 else%first
                      for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=It_left(Delaunay1(f,1),Delaunay1(f,2));
                cf2=It_bottom(Delaunay1(f,1),Delaunay1(f,2));
                cf3=Delaunay1(f,9);
                if cd<cf1 || cd<cf2 ||cd<cf3
                    cc(Delaunay1(f,1),Delaunay1(f,2))=1;
                    
                  
                else
                    
                   cc(xl1+i_-1,yl1+j_-1)=1;
                end
                  
                  
               
                  end
                 end
                           xc=[Delaunay1(f,1)+0.5+v_middle(Delaunay1(f,1),Delaunay1(f,2)),Delaunay1(f,1)+0.5+v_bottom(Delaunay1(f,1),...
                      Delaunay1(f,2)),Delaunay1(f,1)+v_left(Delaunay1(f,1),Delaunay1(f,2))];
                      yc=[Delaunay1(f,2)+0.5+u_middle(Delaunay1(f,1),Delaunay1(f,2)),Delaunay1(f,2)+u_bottom(Delaunay1(f,1),...
                      Delaunay1(f,2)),Delaunay1(f,2)+0.5+u_left(Delaunay1(f,1),Delaunay1(f,2))];
                     [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                     [ii,jj]=find(IN-ON);
                     if size(ii,1)==0
                    
                     else%tird
                      
                      for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=It_left(Delaunay1(f,1),Delaunay1(f,2));
                cf2=It_bottom(Delaunay1(f,1),Delaunay1(f,2));
                cf3=It_middle(Delaunay1(f,1),Delaunay1(f,2));
               if cd<cf1 || cd<cf2 ||cd<cf3
                    cc(Delaunay1(f,1),Delaunay1(f,2))=1;
                    cc( Delaunay1(f,12),Delaunay1(f,13))=1;
                    cc( Delaunay1(f,14),Delaunay1(f,15))=1;
                   
                  
                else
                    cc(xl1+i_-1,yl1+j_-1)=1;
                end
                  
                  
               
                     end
                     end
                            xc=[Delaunay1(f,4),Delaunay1(f,1)+v_left(Delaunay1(f,1),Delaunay1(f,2)),...
                             Delaunay1(f,1)+0.5+v_middle(Delaunay1(f,1),Delaunay1(f,2))];
                             yc=[Delaunay1(f,7),Delaunay1(f,2)+0.5+u_left(Delaunay1(f,1),...
                             Delaunay1(f,2)),Delaunay1(f,2)+0.5+u_middle(Delaunay1(f,1),Delaunay1(f,2))];
                             [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                             [ii,jj]=find(IN-ON);
                             if size(ii,1)==0
                    
                             else%four
                                          
                      for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=It_left(Delaunay1(f,1),Delaunay1(f,2));
               cf2=Delaunay1(f,10);
                cf3=It_middle(Delaunay1(f,1),Delaunay1(f,2));
                 if cd<cf1 || cd<cf2 ||cd<cf3
                        cc( Delaunay1(f,12),Delaunay1(f,13))=1;
                    
                  
                 else
                         cc(xl1+i_-1,yl1+j_-1)=1;
                    
                    
                end
                  
                  
               
                          end
                             end  
                                              xc=[Delaunay1(f,5),Delaunay1(f,1)+0.5+v_bottom(Delaunay1(f,1),Delaunay1(f,2)),...
                                     Delaunay1(f,1)+0.5+v_middle(Delaunay1(f,1),Delaunay1(f,2))];
                                     yc=[Delaunay1(f,8),Delaunay1(f,2)+u_bottom(Delaunay1(f,1),...
                                     Delaunay1(f,2)),Delaunay1(f,2)+0.5+u_middle(Delaunay1(f,1),Delaunay1(f,2))];
                                     [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                                     [ii,jj]=find(IN-ON); 
                                     if size(ii,1)==0
                                     else
                                          for dd=1:size(ii,1)
                             i_=ii(dd);
                               j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=Delaunay1(f,11);
                cf2=It_bottom(Delaunay1(f,1),Delaunay1(f,2));
                cf3=It_middle(Delaunay1(f,1),Delaunay1(f,2));
                if cd<cf1 || cd<cf2 ||cd<cf3
                   cc( Delaunay1(f,14),Delaunay1(f,15))=1;
                    
                  
                else
                   cc(xl1+i_-1,yl1+j_-1)=1;
                    
                end
                  
                  
               
                                           end
                                         
                                     end
                %second
                    else%2
                            xc=[Delaunay1(f,3),Delaunay1(f,1)+v_left(Delaunay1(f,1), Delaunay1(f,2)),Delaunay1(f,1)-0.5+v_middle(Delaunay1(f,1)-1,Delaunay1(f,2))];
                yc=[Delaunay1(f,6),Delaunay1(f,2)+0.5+u_left(Delaunay1(f,1),...
                Delaunay1(f,2)),Delaunay1(f,2)-0.5+u_middle(Delaunay1(f,1)-1,Delaunay1(f,2))];
                [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                [ii,jj]=find(IN-ON);
                 if size(ii,1)==0
           
                 
                 else%first
                      for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=It_left(Delaunay1(f,1),Delaunay1(f,2));
                cf2=It_middle(Delaunay1(f,1)-1,Delaunay1(f,2));
                cf3=Delaunay1(f,9);
                 if cd<cf1 || cd<cf2 ||cd<cf3
                    cc(Delaunay1(f,1),Delaunay1(f,2))=1;
                    
                  
                else
                    
                    cc(xl1+i_-1,yl1+j_-1)=1;
                end
                  
                  
               
                  end
                 end
                            xc=[Delaunay1(f,1)+v_left(Delaunay1(f,1),Delaunay1(f,2)),Delaunay1(f,1)-0.5+v_middle(Delaunay1(f,1)-1,...
                      Delaunay1(f,2)),Delaunay1(f,1)-0.5+v_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1)];
                      yc=[Delaunay1(f,2)+0.5+u_left(Delaunay1(f,1),Delaunay1(f,2)),Delaunay1(f,2)-0.5+u_middle(Delaunay1(f,1)-1,...
                      Delaunay1(f,2)),Delaunay1(f,2)+1+u_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1)];
                     [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                     [ii,jj]=find(IN-ON);
                     if size(ii,1)==0
                   
                     else%tird
                      
                      for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=It_left(Delaunay1(f,1),Delaunay1(f,2));
                cf2=It_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1);
                cf3=It_middle(Delaunay1(f,1)-1,Delaunay1(f,2));
                if cd<cf1 || cd<cf2 ||cd<cf3
                    cc(Delaunay1(f,1),Delaunay1(f,2))=1;
                    cc( Delaunay1(f,12),Delaunay1(f,13))=1;
                    cc( Delaunay1(f,14),Delaunay1(f,15))=1;
                   
                  
                else
                    cc(xl1+i_-1,yl1+j_-1)=1;
                end
                  
                  
               
                     end
                     end
                               xc=[Delaunay1(f,4),Delaunay1(f,1)+v_left(Delaunay1(f,1),Delaunay1(f,2)),...
                             Delaunay1(f,1)-0.5+v_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1)];
                             yc=[Delaunay1(f,7),Delaunay1(f,2)+0.5+u_left(Delaunay1(f,1),...
                             Delaunay1(f,2)),Delaunay1(f,2)+1+u_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1)];
                             [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                             [ii,jj]=find(IN-ON);
                             if size(ii,1)==0
                 
                             else%four
                                          
                      for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=It_left(Delaunay1(f,1),Delaunay1(f,2));
               cf2=Delaunay1(f,10);
                cf3=It_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1);
                 if cd<cf1 || cd<cf2 ||cd<cf3
                     cc( Delaunay1(f,12),Delaunay1(f,13))=1;
                    
                  
                else
                    
                   
                    cc(xl1+i_-1,yl1+j_-1)=1;
                   
                end
                  
                  
               
                          end
                             end
                                                 xc=[Delaunay1(f,5),Delaunay1(f,1)-0.5+v_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1),...
                                     Delaunay1(f,1)-0.5+v_middle(Delaunay1(f,1)-1,Delaunay1(f,2))];
                                     yc=[Delaunay1(f,8),Delaunay1(f,2)+1+u_bottom(Delaunay1(f,1)-1,...
                                     Delaunay1(f,2)+1),Delaunay1(f,2)+0.5+u_middle(Delaunay1(f,1)-1,Delaunay1(f,2))];
                                     [IN,ON]=inpolygon(y2(xl1:xl2,yl1:yl2),x2(xl1:xl2,yl1:yl2),xc,yc);
                                     [ii,jj]=find(IN-ON); 
                                     if size(ii,1)==0
                                     else
                                          for dd=1:size(ii,1)
                             i_=ii(dd);
                               j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
                cf1=Delaunay1(f,11);
                cf2=It_bottom(Delaunay1(f,1)-1,Delaunay1(f,2)+1);
                cf3=It_middle(Delaunay1(f,1)-1,Delaunay1(f,2));
                if cd<cf1 || cd<cf2 ||cd<cf3

                   cc( Delaunay1(f,14),Delaunay1(f,15))=1;
                  
                else
                    
                    
                     cc(xl1+i_-1,yl1+j_-1)=1;
                end
                  
                  
               
                                           end
                                         
                                     end
                %second
                end
            else
               for dd=1:size(ii,1)
                i_=ii(dd);
                j_=jj(dd);
                cd=abs(It(xl1+i_-1,yl1+j_-1));
               if cd<cf1 || cd<cf2 ||cd<cf3
                   cc(Delaunay1(f,1),Delaunay1(f,2))=1;
                    cc( Delaunay1(f,12),Delaunay1(f,13))=1;
                    cc( Delaunay1(f,14),Delaunay1(f,15))=1;
                  
                else
                    
                     cc(xl1+i_-1,yl1+j_-1)=1;
                end
               
            end
            end
           
            
                
            end
            
        end
        figure
 imshow(cc);
 if nargin == 5
count=cc(index).*double(oc(index));
count=length(find(count));
percent =(1- count/countf) *100;   
fprintf('done! missing percent:%.1f \n',percent);

 end


end

