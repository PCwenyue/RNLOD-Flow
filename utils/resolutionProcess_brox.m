function [du, dv] = resolutionProcess_brox(Ikz, Ikx, Iky, omega, uinit, vinit, outer_iter, inner_iter,...
    u_descr,v_descr,ro_descr,delta_inds,lambdaHOG,para,images) 
[ht, wt] = size( Ikz ) ;
if nargin < 8
    %% number of inner iteration
	inner_iter = 500 ;
end
if nargin < 7
    %% number of outer iteration
	outer_iter = 3 ;
end
if nargin < 6
	u = zeros( ht, wt ) ;
	v = zeros( ht, wt ) ;
else
	u = uinit ;
	v = vinit ;
end
if nargin < 4
	omega = 1.9 ;
end
if para.nsigma 
 OCC=para.OCC;
end
du = zeros( ht, wt ) ;
dv = zeros( ht, wt ) ;
tol = 1e-8 * ones( 2 * ht * wt, 1 ) ;
duv = zeros( 2 * ht * wt, 1 ) ;
Psi_MA=zeros( 2 * ht * wt, 1 ) ;
Psi_MB=zeros( 2 * ht * wt, 1 ) ;
Ixx=para.Ixx;Ix2y=para.Ixy;Iyy=para.Iyy;Ixz=para.Ixz;Iyz=para.Iyz;gamma=para.gamma;%% Ixx Ixy Iyy Ixt Iyt gamma
IIxx=Ixx.^2;IIx2y=Ix2y.^2;IIyy=Iyy.^2;
   if size(images,3)==1
              [Ix1 Iy1]= gradient(images(:,:,1));
          else
               [Ix1 Iy1]= gradient(images(:,:,:,1));
                Ix1=mean(Ix1,3);
                Iy1=mean(Iy1,3);
   end
   
   e1=Ix1==0;
   e2=Iy1==0;
   e1=e1|e2;
   Ix1(e1)=1;Iy1(e1)=0;
   
   alg=para.alg;
  
                alpaG=exp(-0.001*(Ix1.^2+Iy1.^2).^0.5);
              
                 IxxG=Ix1./sqrt(Ix1.^2+Iy1.^2);
                 IyyG=Iy1./sqrt(Ix1.^2+Iy1.^2);
                 
                   uv=cat(3,u,v);
                                    F.type=@charbonnier;
                                    F.param= 1.0000e-003;
                                                              
                                    sz        = [size(Ikx,1) size(Ikx,2)];
                                     npixels   = prod(sz);
                                   
for i = 1 : outer_iter
%% HOG
     if size(u_descr,1)~=1
    Psi_descr=get_psi_descr(u(:)+du(:),v(:)+dv(:),u_descr,v_descr,delta_inds);
    Psi_descr=Psi_descr.*ro_descr;
    if para.broxplus2==false
    Psi_MA(1:2:end)=Psi_descr;
    Psi_MA(2:2:end)=Psi_descr;
    Psi_MB(1:2:end)=Psi_descr.*(u_descr-u(:));
    Psi_MB(2:2:end)=Psi_descr.*(v_descr-v(:));
    else
    Psi_MA(1:ht*wt)=Psi_descr;
    Psi_MA(ht*wt+1:end)=Psi_descr;
    Psi_MB(1:ht*wt)=Psi_descr.*(u_descr-u(:));
    Psi_MB(ht*wt+1:end)=Psi_descr.*(v_descr-v(:));    
    end
     end        
   
   Ux=imfilter(u+du,[-1 1 0],'corr','replicate','same');
   Uy=imfilter(u+du,[-1 1 0]','corr','replicate','same');
   Vx=imfilter(v+dv,[-1 1 0],'corr','replicate','same');
   Vy=imfilter(v+dv,[-1 1 0]','corr','replicate','same');
   FU=sparse(npixels, npixels);
   FV=sparse(npixels, npixels);
   rx=deriv_over_x(F,alpaG.*( IxxG.*Ux+IyyG.*Uy).^2+(IyyG.*Ux-IxxG.*Uy).^2,3);
   ry=deriv_over_x(F,alpaG.*( IxxG.*Vx+IyyG.*Vy).^2+(IyyG.*Vx-IxxG.*Vy).^2,3);
for iG=1:3
   if iG==1
     FMi = make_convn_mat([1 -1], sz, 'valid', 'sameswap');
     Fi=FMi';
   temp=(alpaG.*IxxG.^2+IyyG.^2).*rx;
   
   Fu        = Fi*spdiags(temp(:), 0, npixels, npixels)*FMi;
    temp=(alpaG.*IxxG.^2+IyyG.^2).*ry;
   Fv        = Fi*spdiags(temp(:), 0, npixels, npixels)*FMi;
       FU=FU-Fu;
       FV=FV-Fv;
   elseif iG==2
          FMi = make_convn_mat([1; -1], sz, 'valid', 'sameswap');
     Fi=FMi';
   temp=(alpaG.*IyyG.^2+IxxG.^2).*rx;
   Fu        =Fi*spdiags(temp(:), 0, npixels, npixels)*FMi;
   temp=(alpaG.*IyyG.^2+IxxG.^2).*ry;
   Fv        =Fi*spdiags(temp(:), 0, npixels, npixels)*FMi;
       FU=FU-Fu;
       FV=FV-Fv;
   else
       FMi = make_imfilter_mat([1 0 -1], sz, 'replicate', 'same');
       FMI= make_imfilter_mat([1 0 -1]', sz, 'replicate', 'same');
        temp=1/2*(alpaG-ones(size(IxxG))).*IxxG.*IyyG.*rx;
        Fu        = FMI*spdiags(temp(:), 0, npixels, npixels)*FMi;
          temp=1/2*(alpaG-ones(size(IxxG))).*IxxG.*IyyG.*ry;
        Fv        = FMI*spdiags(temp(:), 0, npixels, npixels)*FMi;
         FU=FU+Fu;
       FV=FV+Fv;
   end
end  
 M = [FU, sparse(npixels, npixels);
      sparse(npixels, npixels),FV];                           
                        
                          Ix2 = Ikx.^2;
                          Iy2 = Iky.^2;
                          Ixy = Ikx.*Iky;
                          Itx = Ikz.*Ikx;
                          Ity = Ikz.*Iky;                

                          It = (Ikz + Ikx.*repmat(du, [1 1 size(Ikz,3)]) ...
                          + Iky.*repmat(dv, [1 1 size(Ikz,3)])).^2+gamma*((Ixz+Ixx.*du+Ix2y.*dv).^2+(Iyz+Ix2y.*du+Iyy.*dv).^2);
                      if para.nsigma         
                            pp_d =OCC(:).*deriv_over_x(F, It(:),3);
                      else
                            pp_d =deriv_over_x(F, It(:),3);
                      end

                           tmp = pp_d.*(Ix2(:)+gamma*(IIxx(:)+IIx2y(:)));
                          duu = spdiags(tmp, 0, npixels, npixels);

                          tmp = pp_d.*(Iy2(:)+gamma*(IIyy(:)+IIx2y(:)));
                          dvv = spdiags(tmp, 0, npixels, npixels);

                          tmp = pp_d.*(Ixy(:)+gamma*(Ixx(:).*Ix2y(:)+Iyy(:).*Ix2y(:)));
                          dduv = spdiags(tmp, 0, npixels, npixels);
                         
                          A = [duu dduv; dduv dvv] - M*alg;
                         
                           b = alg*M*uv(:) - [pp_d.*(Itx(:)+gamma*(Ixx(:).*Ixz(:)+Ix2y(:).*Iyz(:))); pp_d.*(Ity(:)+gamma*(Ix2y(:).*Ixz(:)+Iyy(:).*Iyz(:)))];
                           
                          if lambdaHOG~=0 
                          A= A+spdiags(lambdaHOG*Psi_MA, 0, size(A,1), size(A,2));
                          b= b + lambdaHOG*Psi_MB; 
                          end
 	[duv, err, it, flag] = sor( A, duv, b, omega, inner_iter, tol ) ;%% SOR
                    duv(duv > 1)  = 1;
                    duv(duv < -1) = -1;            
                       du(:)= duv( 1:npixels )  ;
                       dv(:)= duv( npixels+1:end )  ;
end
 
        
