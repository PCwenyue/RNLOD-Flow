%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%   This software is provided for research purposes only. Any commercial use is prohibited. 
%   Any scientific work that makes use of our code should cite our published paper: 
%   Zhang Congxuan, Chen Zhen, Wang Mingrun, Li Ming, Jaing Shaofeng, "Robust Non-Local TV-L1 Optical Flow Estimation with Occlusion Detection", IEEE Transactions on Image Processing, 2017, DOI: 10.1109/TIP.2017.2712279
%   This version of the code is updated and not fixed the optimized parameters, you can adjust the parameters in the demo file to reach the better performance.  
%   Date: May 5, 2017, Editor: Zhang Congxuan, Nanchang Hangkong University
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
addpath('utils');
%%   ===============================================
 ima1={'Dimetrodon','Grove2','Grove3','Hydrangea','Urban2','Urban3','Venus','Rubberwhale'};
%%  ================================================
para.gamma=0.5;
para.lambdaHOG=0;  % scale of HOG, 0 denotes No, 1 denotes Yes
para.blend=0.95;   % parameter of structure_texture_decomposition
para.num_levels=30;% Layers of pyramid
para.smooth_sigma=0.5; % standard deviation of Gaussian filter 
para.space=0.9;    % scale of the pyramid down-sampling 
para.h=[-1/2 0 1/2]; % derivation template
para.alg=2;        % scale of the smoothing term
para.OCC=[];       % initial occlusion region, set as []
para.nsigma=false;  
para.fullversion=true;
para.alphaG=0.5;

%% =====================================================
%%% middlebury_other-data_training_sequence
for  i=2
 
filename=ima1{i};
imagename='frame10';
imagename2='frame11';
img1=imread(['train_data\other-color-data\',filename,'\',imagename,'.png']);
img2=imread(['train_data\other-color-data\',filename,'\',imagename2,'.png']);
imagename='flow10';   
flow1 = readFlowFile(['train_data\other-gt-flow\',filename,'\',imagename,'.flo']);
[u v Rtime]=optic_flow_compute(img1, img2,para);                        % main function
display11(u,v,flow1(:,:,1),flow1(:,:,2),filename,Rtime);
occusiondetection( cat(3,u,v),cat(4,double(img1),double(img2)));
% writeFlowFile(cat(3,u,v), filename);
end

















