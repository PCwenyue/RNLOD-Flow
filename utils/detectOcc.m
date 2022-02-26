function [  div It1 OCC] = detectOcc(uv, images)
%% compute Edata£¨It£©
It1 = partial_deriv(images, uv);
%% divergence of optical flow
 div = divergence(uv(:,:,1), uv(:,:,2));
 if nargout==3
       bfhsz=10; 
It=abs(It1);       
%% Initialization of occlusion detection
H=size(uv,1);W=size(uv,2);
                    OCC=false(H,W);
                    fl=ones(H,W);
                    fl(end,:)=0;
                    fl(:,end)=0;
                    index=find(fl==0);
                    [C R]=meshgrid(1:W,1:H);
                    U=uv(:,:,1)+C;
                    V=uv(:,:,2)+R;
                    %% Define the region size of detection
                    e1 = edge(uv(:,:,1), 'sobel');
                    e2 = edge(uv(:,:,2), 'sobel');
                    e  = e1|e2;
                    mask = imdilate(e, ones([10 10]) );
                       
                  
                    count=find(mask);

                   
                    [c r] = meshgrid(-bfhsz:bfhsz, -bfhsz:bfhsz);
                   
                         i=1;
                    while i<length(count)+1
                        %% Prevent the search templates over the image boundary
                            if any(count(i)==index)
                                i=i+1;
                                continue;
                            end
                           if rem(count(i),H)-1 < bfhsz
                               r1=r(bfhsz+2-rem(count(i),H):end,:);
                                c1=c(bfhsz+2-rem(count(i),H):end,:);
                           elseif rem(count(i),H)+bfhsz>H
                               r1=r(1:end-rem(count(i),H)-bfhsz+H,:);
                               c1=c(1:end-rem(count(i),H)-bfhsz+H,:);
                           else
                               r1=r;
                               c1=c;
                           end
                           [k l]=ind2sub([H,W],count(i));
                           if l-1 < bfhsz
                               r1=r1(:,bfhsz-l+2:end);
                               c1=c1(:,bfhsz-l+2:end);
                           elseif W-l < bfhsz
                               r1=r1(:,1:end-bfhsz+W-l);
                               c1=c1(:,1:end-bfhsz+W-l);
                           end
                           %% Calculate coordinates of each pixel in the second frame
                           nindx = r1 + c1*H;
                            x1=U(count(i)); 
                            x2=U(count(i)+H); 
                            x3=U(count(i)+1);
                            y1=V(count(i));
                            y2=V(count(i)+H); 
                            y3=V(count(i)+1);
                            %% Determine whether the pixel is embedded in the triangle
                            [IN ON]=inpolygon(U(nindx(:)+count(i)),V(nindx(:)+count(i)),[x1 x2 x3]',[y1 y2 y3]');
                          
                            j =find(IN-ON);
                     
                            for ii=1:length(j)
                                %% distinguish whether the triangle or the embedded pixel is occluded
                                  if It(nindx(j(ii))+count(i))>(It(count(i))+It(count(i)+H)+It(count(i)+1))/3
                                      
                                OCC(nindx(j(ii))+count(i))=true;
                                  else
                                      
                                       OCC(count(i))=true;   
                                       OCC(count(i)+H)=true;
                                       OCC(count(i)+1)=true;

                                       
                                 end
     
                            end
                           
                            i=i+1;
                           
                    end
 end

