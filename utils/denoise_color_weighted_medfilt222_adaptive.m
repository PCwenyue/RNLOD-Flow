function  uvo = denoise_color_weighted_medfilt222_adaptive(uv, im,div,it,OCC,fullversion,ap)
if nargin<6
    OCC=false(size(it));   
end

          
 %%  div>0
level = graythresh(div(div>0));
 %%  div<0
level2 = graythresh(abs(div(div<0)));
nIRLS=5; 
sz = size(im);
sz = sz(1:2);
f = fspecial('gaussian',5, 0.5) ;  
    bfhsz =10; 
    mfhsz=bfhsz;
    alpha=0.05;   
    TR = 15;
    Tau = 0.7;
    Tau_uv=1;
    tau=2;
    M = sqrt(uv(:,:,1).^2 + uv(:,:,2).^2);
    RMS_M =mean(M(:));
  
        %% determine the weight of median filter
        sigma_c=12;
        sigma_i =3; 
        sigma_d =0.3; 
        sigma_x =7; 
        %%%%%
         sigma_c_occ=3;
        sigma_i_occ =3; 
        sigma_d_occ =0.3; 
        sigma_x_occ =7; 
        %%%%%
        sigma_c_0=5;
         sigma_i_0=3; 
        sigma_d_0 =0.3; 
        sigma_x_0 =7; 
        %%%%%
        sigma_c_3=3;
        sigma_d_3=0.3;
        sigma_x_3 =7; 



 divn=div>level; 
 divp=div<-level2; 
 uvo=uv;
 ind=find(divp);
 ind0=find(divn);
 indocc=find(~(divp|divn)&OCC);

div(div>0)=0;


if fullversion
   mask = ones(size(div)); 
else
e1 = edge(uv(:,:,1), 'sobel');
e2 = edge(uv(:,:,2), 'sobel');
e  = e1|e2;
mask =imdilate(e, ones([5 5]) );
end
    


[indx_row, indx_col] = find(mask );

im = imfilter(im, f, 'corr', 'symmetric', 'same');  
it = imfilter(it, f, 'corr', 'symmetric', 'same'); 

  levelit = graythresh(abs(it));
  mask_it=abs(it)>levelit;

pad_it=padarray(it, bfhsz*[1 1], 'symmetric', 'both'); 
pad_mask_it=padarray(mask_it, bfhsz*[1 1], 'symmetric', 'both'); 
pad_div= padarray(div, bfhsz*[1 1], 'symmetric', 'both');     
pad_u  = padarray(uvo(:,:,1), bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(uvo(:,:,2), bfhsz*[1 1], 'symmetric', 'both');        
pad_im = padarray(im, bfhsz*[1 1], 'symmetric', 'both');    
pad_divp= padarray(divp, bfhsz*[1 1], 'symmetric', 'both');


    rho.type='generalized_charbonnier';
    rho.param=[1e-3, 0.45];


[H W] = size(pad_u);

% Divide into several groups for memory reasons ~70,000 causes out of memory
Indx_Row = indx_row;
Indx_Col = indx_col;
N        = length(Indx_Row); % number of elements to process
n        = 4e4;              % number of elements per batch
nB       = ceil(N/n);

for ib = 1:nB;
   
    istart = (ib-1)*n + 1;
    iend   = min(ib*n, N);
    ind1= sub2ind(sz, Indx_Row(istart:iend), Indx_Col(istart:iend));
    ind3=ismember(ind1,ind);
    ind0=ismember(ind1,ind0);
     indocc=ismember(ind1,indocc);
    ind4=~(ind3+ind0+indocc);
    indx_row = Indx_Row(istart:iend);
    indx_col = Indx_Col(istart:iend);    

    [C R] = meshgrid(-bfhsz:bfhsz, -bfhsz:bfhsz);
    
    nindx = R + C*H;    

    cindx = indx_row +bfhsz  + (indx_col+bfhsz-1)*H;
    %%  
    cin=cindx(ind4);
    pad_indx = repmat(nindx(:), [1 length(cin)]) + ...
               repmat(cin(:)', [(bfhsz*2+1)^2, 1] );
    
    %%  
     tmp = exp(- (C.^2 + R.^2) /2/sigma_x^2 );
 
%      e1=tmp(:);
     weights = repmat(tmp(:), [1 length(cin)]);    
  
%       sd = repmat(tmp1(:), [1 length(indx_row)]);    
%       weights1=sd;
    % %Uncomment below: no spatial weight for test
    % weights = ones(size(weights));    
  
     tmp_w = zeros(size(weights));
   
    for i = 1:size(pad_im,3)
        tmp = pad_im(:,:,i);
        tmp_w = tmp_w + (tmp(pad_indx) - repmat(tmp(cin(:))', [(bfhsz*2+1)^2, 1])).^2;
    end;    
    tmp_w = tmp_w/size(pad_im,3);

     weights = weights.* exp(-tmp_w/2/sigma_c^2);
    weights=weights.*exp(-pad_div(pad_indx).^2/2/sigma_d^2).*exp(-pad_it(pad_indx).^2/2/sigma_i^2);
   
    % Occlusion weight    
     
%      weights(SE,:)=[];
    
%       weights(:,ind3)=ext(:,ind3);
  
   
  
    weights = weights./repmat(sum(weights, 1), [(2*bfhsz+1)^2, 1]);
    
    neighbors_u = pad_u(pad_indx);
    neighbors_v = pad_v(pad_indx);

      uo   = uvo(:,:,1);
      u    = weighted_median(weights, neighbors_u);

    uo(ind1(ind4))=u;
    vo   = uvo(:,:,2);
    v    = weighted_median(weights, neighbors_v);

vo(ind1(ind4))=v;
%% uplate uv
pad_u  = padarray(uo, bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(vo, bfhsz*[1 1], 'symmetric', 'both');   
%% occlusion 
cin=cindx(indocc);
if ~isempty(cin)
    pad_indx = repmat(nindx(:), [1 length(cin)]) + ...
               repmat(cin(:)', [(bfhsz*2+1)^2, 1] );
         
        
           
           temp_u=zeros(size(cin));
           temp_v=zeros(size(cin));
            
        
           
           %%no it
            
        
    % spatial weight
    te=(C.^2 + R.^2) /2/sigma_x_3^2 ;
     tmp = exp(- te);
 
%      e1=tmp(:);
     weights = repmat(tmp(:), [1 length(cin)]); 
     
       tmp_w = zeros(size(weights));
       tmp_w1 = tmp_w;
    for i = 1:size(pad_im,3)
        tmp = pad_im(:,:,i);
        tmpw=tmp(pad_indx) - repmat(tmp(cin(:))', [(bfhsz*2+1)^2, 1]);
        tmp_w = tmp_w + (tmpw).^2;
        tmp_w1 = tmp_w1 + abs(tmpw);
    end;    
    tmp_w = tmp_w/size(pad_im,3);
    tmp_w1 = tmp_w/size(pad_im,3);
    %%%%%%%%%s
    
%% Motion Distance Difference
    % Motionu = (pad_u(pad_indx) - repmat(pad_u(cindx(:))', [(mfhsz*2+1)^2, 1]));     % abs(.^1); .^2
    % Motionv = (pad_v(pad_indx) - repmat(pad_v(cindx(:))', [(mfhsz*2+1)^2, 1]));         
    Padu = pad_u(pad_indx);
    PadReu = repmat(pad_u(cin(:))', [(mfhsz*2+1)^2, 1]);
    Padv = pad_v(pad_indx);
    PadRev = repmat(pad_v(cin(:))', [(mfhsz*2+1)^2, 1]);
    Motionu = abs(Padu-PadReu);
    Motionv = abs(Padv-PadRev);
    
%% Find the ponits which with small motion difference(u,v) while large intensity difference
    Meanu = mean((Motionu));                      % mean, median
    Meanv = mean((Motionv));
    Meani = mean(tmp_w1);  
    Meanu = repmat(Meanu, [(mfhsz*2+1)^2, 1]);
    Meanv = repmat(Meanv, [(mfhsz*2+1)^2, 1]);
    Meani = repmat(Meani, [(mfhsz*2+1)^2, 1]);
     
    Labeliu = ((Motionu) > Meanu*Tau_uv)&(tmp_w1 > Tau*Meani);          % u diff very small(similar motion patch)& color differ very large;1/2 
    Labeliv = ((Motionv) > Meanv*Tau_uv)&(tmp_w1 > Tau*Meani);
%% Motion Angle Difference (Large displacements)
 
          if RMS_M>=TR  
        Anguv = atand(Padu./(Padv+0.001));
        AngRuv = atand(PadReu./(PadRev+0.001));   
        AngS = Anguv.*AngRuv;
        AngD = abs(Anguv - AngRuv);
        Labeliu = (Labeliu)&(AngS>=0)&(AngD<=TR);         
        Labeliv = (Labeliv)&(AngS>=0)&(AngD<=TR);
          end
           Labeliuv = Labeliu|Labeliv;    
  
           sumn=sum( Labeliuv ,1);
%            sumn(sumn==0)=1;
                 
            tmpp=pad_mask_it(pad_indx);
            tmpp1=pad_divp(pad_indx);
            fa=~(tmpp|Labeliuv|tmpp1);
            sump=sum(fa,1);
           
           ind_noit=sump<size(pad_indx,1)*alpha;    
           ind_noG=~ind_noit&(sumn==0);
           ind_it=xor(~ind_noit,ind_noG);
    tew=tmp_w/2/sigma_c_3^2;
    tw=exp(-tew);

  %% compute factor I
 
    weights(:,ind_noit) = weights(:,ind_noit).*tw(:,ind_noit).*exp(-pad_div(pad_indx(:,ind_noit)).^2./2/sigma_d_3^2) ;
    weights(:,ind_noG) = weights(:,ind_noG).*tw(:,ind_noG).*exp(-pad_div(pad_indx(:,ind_noG)).^2./2/sigma_d_occ^2).*exp(-pad_it(pad_indx(:,ind_noG)).^2./2/sigma_i_occ^2) ;
            if sum(ind_it)~=0
    exn=sum(repmat(te(:), [1 sum(ind_it)]).*Labeliuv(:,ind_it) ,1)./sumn(ind_it);
      
        expp=sum(repmat(te(:), [1 sum(ind_it)]).*fa(:,ind_it),1)./sump(ind_it);
   
    %%
    cn=sum(tew(:,ind_it).*Labeliuv(:,ind_it))./sumn(ind_it);
  
    I=pad_it(cin(ind_it)').^2./(ap*cn+ap*exn+(1-ap)*expp);
    I=repmat(I, [(2*bfhsz+1)^2, 1]);
    I(I==0)=0.001;
   weights(:,ind_it) = weights(:,ind_it).*tw(:,ind_it).*exp(-pad_div(pad_indx(:,ind_it)).^2./I).*exp(-pad_it(pad_indx(:,ind_it)).^2./I);
            end
    % Normalize
    weights = weights./repmat(sum(weights, 1), [(2*bfhsz+1)^2, 1]);
     weights(isnan(weights))=1/size(pad_indx,1);
    neighbors_u = pad_u(pad_indx);
    neighbors_v = pad_v(pad_indx);
    weightsColor = weights; 
    u = sum(weights.*neighbors_u);
    v = sum(weights.*neighbors_v);

    for iter=1:nIRLS
        % iterated least square: change rho to a weighted L2
        tmp_w = pad_u(pad_indx) - repmat(u, [(bfhsz*2+1)^2, 1]);
        tmp_w = deriv_over_x(rho, tmp_w);
        weights = weightsColor.*tmp_w;
        
        tmp_w = pad_v(pad_indx) - repmat(v, [(bfhsz*2+1)^2, 1]);
        tmp_w = deriv_over_x(rho, tmp_w);
        weights = weights.*tmp_w;
        
        % Normalize
        weights = weights./repmat(sum(weights, 1), [(bfhsz*2+1)^2, 1]);
        u = sum(weights.*neighbors_u);
        v = sum(weights.*neighbors_v);
    end;

  

   
   
   
uo(ind1(indocc))=u;
vo(ind1(indocc))=v;

%% uplate uv
pad_u  = padarray(uo, bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(vo, bfhsz*[1 1], 'symmetric', 'both');   
end
%%%%%%%%%%%%%%%%%%%%%%
  cin=cindx(ind0);
  if ~isempty(cin)
    pad_indx = repmat(nindx(:), [1 length(cin)]) + ...
               repmat(cin(:)', [(bfhsz*2+1)^2, 1] );
    
    % spatial weight
     tmp = exp(- (C.^2 + R.^2) /2/sigma_x_0^2 );
 
%      e1=tmp(:);
     weights = repmat(tmp(:), [1 length(cin)]);    
  
%       sd = repmat(tmp1(:), [1 length(indx_row)]);    
%       weights1=sd;
    % %Uncomment below: no spatial weight for test
    % weights = ones(size(weights));    
    % Intensity weight
     tmp_w = zeros(size(weights));
   
    for i = 1:size(pad_im,3)
        tmp = pad_im(:,:,i);
        tmp_w = tmp_w + (tmp(pad_indx) - repmat(tmp(cin(:))', [(bfhsz*2+1)^2, 1])).^2;
    end;    
    tmp_w = tmp_w/size(pad_im,3);

     weights = weights.* exp(-tmp_w/2/sigma_c_0^2);
    weights=weights.*exp(-pad_div(pad_indx).^2/2/sigma_d_0^2).*exp(-pad_it(pad_indx).^2/2/sigma_i_0^2);
   
    % Occlusion weight    
     
%      weights(SE,:)=[];
    
%       weights(:,ind3)=ext(:,ind3);
  
   
    % Normalize
    weights = weights./repmat(sum(weights, 1), [(2*bfhsz+1)^2, 1]);
    
    neighbors_u = pad_u(pad_indx);
    neighbors_v = pad_v(pad_indx);

    weightsColor = weights; 
    u = sum(weights.*neighbors_u);
    v = sum(weights.*neighbors_v);

    for iter=1:nIRLS
        % iterated least square: change rho to a weighted L2
        tmp_w = pad_u(pad_indx) - repmat(u, [(bfhsz*2+1)^2, 1]);
        tmp_w = deriv_over_x(rho, tmp_w);
        weights = weightsColor.*tmp_w;
        
        tmp_w = pad_v(pad_indx) - repmat(v, [(bfhsz*2+1)^2, 1]);
        tmp_w = deriv_over_x(rho, tmp_w);
        weights = weights.*tmp_w;
        
        % Normalize
        weights = weights./repmat(sum(weights, 1), [(bfhsz*2+1)^2, 1]);
        u = sum(weights.*neighbors_u);
        v = sum(weights.*neighbors_v);
    end;

   
uo(ind1(ind0))=u;
vo(ind1(ind0))=v;
pad_u  = padarray(uo, bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(vo, bfhsz*[1 1], 'symmetric', 'both');   
  end
%%%%%%%%%%%%%%%%%%%%%%
 cin=cindx(ind3);
 if ~isempty(cin)
    pad_indx = repmat(nindx(:), [1 length(cin)]) + ...
               repmat(cin(:)', [(bfhsz*2+1)^2, 1] );
      
           
           %%no it
            
       
    % spatial weight
    te=(C.^2 + R.^2) /2/sigma_x_3^2 ;
     tmp = exp(- te);
 
%      e1=tmp(:);
     weights = repmat(tmp(:), [1 length(cin)]);    
              tmp_w = zeros(size(weights));
    tmp_w1= tmp_w;
    for i = 1:size(pad_im,3)
        tmp = pad_im(:,:,i);
        tmp1=tmp(pad_indx) - repmat(tmp(cin(:))', [(bfhsz*2+1)^2, 1]);
        tmp_w = tmp_w + (tmp1).^2;
         tmp_w1 = tmp_w1 + abs(tmp1);
    end;    
    tmp_w = tmp_w/size(pad_im,3);
     tmp_w1 = tmp_w1/size(pad_im,3);
    %%%%%%%%%s
    
%% Motion Distance Difference
    % Motionu = (pad_u(pad_indx) - repmat(pad_u(cindx(:))', [(mfhsz*2+1)^2, 1]));     % abs(.^1); .^2
    % Motionv = (pad_v(pad_indx) - repmat(pad_v(cindx(:))', [(mfhsz*2+1)^2, 1]));         
    Padu = pad_u(pad_indx);
    PadReu = repmat(pad_u(cin(:))', [(mfhsz*2+1)^2, 1]);
    Padv = pad_v(pad_indx);
    PadRev = repmat(pad_v(cin(:))', [(mfhsz*2+1)^2, 1]);
    Motionu = abs(Padu-PadReu);
    Motionv = abs(Padv-PadRev);

%% Find the ponits which with small motion difference(u,v) while large intensity difference
    Meanu = mean((Motionu));                      % mean, median
    Meanv = mean((Motionv));
    Meani = mean(tmp_w1);   
    Meanu = repmat(Meanu, [(mfhsz*2+1)^2, 1]);
    Meanv = repmat(Meanv, [(mfhsz*2+1)^2, 1]);
    Meani = repmat(Meani, [(mfhsz*2+1)^2, 1]);
     
    Labeliu = ((Motionu) < Meanu/Tau)&(tmp_w1 > Tau*Meani)&((Motionu) < tau);          % u diff very small(similar motion patch)& color differ very large;1/2 
    Labeliv = ((Motionv) < Meanv/Tau)&(tmp_w1 > Tau*Meani)&((Motionv) < tau);
%% Motion Angle Difference (Large displacements)

       if RMS_M>=TR  
        Anguv = atand(Padu./(Padv+0.001));
        AngRuv = atand(PadReu./(PadRev+0.001));   
        AngS = Anguv.*AngRuv;
        AngD = abs(Anguv - AngRuv);
        Labeliu = (Labeliu)&(AngS>=0)&(AngD<=TR);         
        Labeliv = (Labeliv)&(AngS>=0)&(AngD<=TR);
    end


    Labeliuv = Labeliu|Labeliv;    
    
       
             sumn=sum( Labeliuv ,1);
%            sumn(sumn==0)=1;
                 
            tmpp=pad_mask_it(pad_indx);
            tmpp1=pad_divp(pad_indx);
            fa=~(tmpp|Labeliuv|tmpp1);
            sump=sum(fa,1);
           
           ind_noit=sump<size(pad_indx,1)*alpha;    
           ind_noG=~ind_noit&(sumn==0);
           ind_it=xor(~ind_noit,ind_noG);
    tew=tmp_w/2/sigma_c_3^2;
    tw=exp(-tew);

  %% 
 
    weights(:,ind_noit) = weights(:,ind_noit).*tw(:,ind_noit).*exp(-pad_div(pad_indx(:,ind_noit)).^2./2/sigma_d_3^2) ;
    weights(:,ind_noG) = weights(:,ind_noG).*tw(:,ind_noG).*exp(-pad_div(pad_indx(:,ind_noG)).^2./2/sigma_d_occ^2).*exp(-pad_it(pad_indx(:,ind_noG)).^2./2/sigma_i_occ^2) ;
            if sum(ind_it)~=0
    exn=sum(repmat(te(:), [1 sum(ind_it)]).*Labeliuv(:,ind_it) ,1)./sumn(ind_it);
            
      
        expp=sum(repmat(te(:), [1 sum(ind_it)]).*fa(:,ind_it),1)./sump(ind_it);
   
    %% 
    cn=sum(tew(:,ind_it).*Labeliuv(:,ind_it))./sumn(ind_it);
  
    I=pad_it(cin(ind_it)').^2./(ap*cn+ap*exn+(1-ap)*expp);
    I=repmat(I, [(2*bfhsz+1)^2, 1]);
    I(I==0)=0.001;
   weights(:,ind_it) = weights(:,ind_it).*tw(:,ind_it).*exp(-pad_div(pad_indx(:,ind_it)).^2./I).*exp(-pad_it(pad_indx(:,ind_it)).^2./I);
            end
    % Normalize
    weights = weights./repmat(sum(weights, 1), [(2*bfhsz+1)^2, 1]);
     weights(isnan(weights))=1/size(pad_indx,1);
    neighbors_u = pad_u(pad_indx);
    neighbors_v = pad_v(pad_indx);
    weightsColor = weights; 
    u = sum(weights.*neighbors_u);
    v = sum(weights.*neighbors_v);

    for iter=1:nIRLS
        % iterated least square: change rho to a weighted L2
        tmp_w = pad_u(pad_indx) - repmat(u, [(bfhsz*2+1)^2, 1]);
        tmp_w = deriv_over_x(rho, tmp_w);
        weights = weightsColor.*tmp_w;
        
        tmp_w = pad_v(pad_indx) - repmat(v, [(bfhsz*2+1)^2, 1]);
        tmp_w = deriv_over_x(rho, tmp_w);
        weights = weights.*tmp_w;
        
        % Normalize
        weights = weights./repmat(sum(weights, 1), [(bfhsz*2+1)^2, 1]);
        u = sum(weights.*neighbors_u);
        v = sum(weights.*neighbors_v);
    end;

  

   
   
   
uo(ind1(ind3))=u;
vo(ind1(ind3))=v;
 end


 

    uvo = cat(3, uo, vo);

end;

