function display1( u,v,tu,tv ,filename,mask)
flow=cat(3,u,v);
 imflow = flowToColor(flow);                                                

 UNKNOWN_FLOW_THRESH = 1e9;
 tu (tu>UNKNOWN_FLOW_THRESH) = NaN;
 tv (tv>UNKNOWN_FLOW_THRESH) = NaN;
 if sum(~isnan(tu(:))) > 1
 
     [aae stdae aepe] = flowAngErr(tu,tv, u, v, 0);
     fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
      [oaae ostdae oaepe] = flowAngErr1(tu,tv, u, v, 0,mask);
     fprintf('\nAAE %3.3f average EPE %3.3f \n', oaae, oaepe);
     
 end;
 straae=num2str(round(aae*1000)/1000);
 strepe=num2str(round(aepe*1000)/1000);
 ostraae=num2str(round(oaae*1000)/1000);
 ostrepe=num2str(round(oaepe*1000)/1000);
 figure
 subplot(1,1,1);imshow(imflow); title([filename,' aae:',straae,' aepe:',strepe,' oaepe:',ostrepe,' oaae:',ostraae]);


end

