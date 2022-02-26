function display11( u,v,tu,tv ,filename,Rtime)
 flow=cat(3,u,v);
 imflow = flowToColor(flow);                                                

 UNKNOWN_FLOW_THRESH = 1e9;
 tu (tu>UNKNOWN_FLOW_THRESH) = NaN;
 tv (tv>UNKNOWN_FLOW_THRESH) = NaN;
 if sum(~isnan(tu(:))) > 1
     [aae, staae, aepe]= flowAngErr(tu,tv, u, v, 0);
%      [aae stdae aepe] = flowAngErr(tu,tv, u, v, 0);
     fprintf('\naae %3.3f average aepe %3.3f \n',aae, aepe);   
 end;
 straae=num2str(round(aae*1000)/1000);
 straepe=num2str(round(aepe*1000)/1000);
 strtime=num2str(round(Rtime*1000)/1000);
 figure
 subplot(1,1,1);imshow(imflow); title([filename,' aae:',straae,' aepe:',straepe,'Runtime:',strtime]);
% figure
%  subplot(1,1,1);imshow(imflow); title([filename,' aae:',aae,' aepe:',aepe]);
% % 
end

