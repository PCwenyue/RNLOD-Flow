function [ uv DOCC] = getMedianuv( uv , img_color ,images,para)

                      [  div it OCC] = detectOcc(uv, images);  % detect the occlusion region


    uv=denoise_color_weighted_medfilt222_adaptive(uv,  img_color,div,it, OCC,para.fullversion,para.alphaG);% median filter

if nargout==2
    DOCC=OCC;
end