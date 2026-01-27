%% very crude translation of 1D heminode address to 2D. Try not to think too much about it :P
% assumes 8x8 hem array
% function [x,y] = linto2D(h)
% 
%     if h >= 1 && h <= 8
%        x = h; 
%        y = 1;
%     elseif h >= 9 && h <= 16
%        x = h - 8;
%        y = 2;
%     elseif h >= 17 && h <= 24
%        x = h - 16;
%        y = 3;
%     elseif h >= 25 && h <= 32
%        x = h - 24;
%        y = 4;
%     elseif h >= 33 && h <= 40
%        x = h - 32;
%        y = 5;
%     elseif h >= 41 && h <= 48
%        x = h - 40;
%        y = 6;
%     elseif h >= 49 && h <= 56
%        x = h - 48;
%        y = 7;
%     elseif h >= 57 && h <= 64
%        x = h - 56;
%        y = 8;       
%     end      
%     
% end


function [x,y] = linto2D(h,mr_cnt)

mr_cnt_oneside = sqrt(mr_cnt);
mr_floor = floor((h-1)/mr_cnt_oneside);

    %if h >= 1 && h <= mr_cnt_oneside
       x = h - mr_floor*mr_cnt_oneside; 
       y = mr_floor+1;
%     elseif h >= 9 && h <= 16
%        x = h - 8;
%        y = 2;
%     elseif h >= 17 && h <= 24
%        x = h - 16;
%        y = 3;
%     elseif h >= 25 && h <= 32
%        x = h - 24;
%        y = 4;
%     elseif h >= 33 && h <= 40
%        x = h - 32;
%        y = 5;
%     elseif h >= 41 && h <= 48
%        x = h - 40;
%        y = 6;
%     elseif h >= 49 && h <= 56
%        x = h - 48;
%        y = 7;
%     elseif h >= 57 && h <= 64
%        x = h - 56;
%        y = 8;       
%     end      
    
end