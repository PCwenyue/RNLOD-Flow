function [ occ] = detecOCC( uv, images )
bfhsz=10;
u=uv(:,:,1);
v=uv(:,:,2);
sigma_d = 0.3; 
sigma_i =0.5;  
                    H=size(uv,1);W=size(uv,2);
                    OCC=ones(H,W);
                    fl=ones(H,W);
                    fl(end,:)=0;
                    fl(:,end)=0;
                    [C R]=meshgrid(1:W,1:H);
                    U=u+C;
                    V=v+R;
                    pad_u  = padarray(u, bfhsz*[1 1], 'replicate', 'both');        
                    pad_v  = padarray(v, bfhsz*[1 1], 'replicate', 'both'); 
                    pad_U  = padarray(U, bfhsz*[1 1], 'replicate', 'both');        
                    pad_V  = padarray(V, bfhsz*[1 1], 'replicate', 'both'); 
                    e1 = edge(uv(:,:,1), 'sobel');
                    e2 = edge(uv(:,:,2), 'sobel');
                    e  = e1|e2;
                    mask = imdilate(e, ones([5 5]) );
                     mask(fl==0)=0;
                    It = partial_deriv(images, uv);
                    [indx_row, indx_col] =find(mask);
                    count=find(mask);
                    [c r] = meshgrid(-bfhsz:bfhsz, -bfhsz:bfhsz);
                    nindx = r + c*H;
                    cindx = indx_row +bfhsz  + (indx_col+bfhsz-1)*H;
                    cindx=repmat(cindx(:)', [(bfhsz*2+1)^2, 1] );
                    pad_indx = repmat(nindx(:), [length(indx_row) 1 ]) + cindx(:);
                    pad_indx=unique( pad_indx);
                    x(1:4:end)=U(count);
                    x(2:4:end)=U(count+H);
                    x(3:4:end)=U(count+1);
                    x(4:4:end)=U(count);
                    y(1:4:end)=V(count);
                    y(2:4:end)=V(count+H);
                    y(3:4:end)=V(count+1);
                    y(4:4:end)=V(count);
                    [IN ON]=inpolygon(pad_U(pad_indx),pad_V(pad_indx),x,y);%内接多边形
                    IN=IN-ON;
                    

              It=It.*OCC;


            div = divergence(uv(:,:,1), uv(:,:,2));
            div(div>0) =0;
            occ = exp(-div.^2/2/sigma_d^2).*exp(-It.^2/2/ sigma_i^2);


end

