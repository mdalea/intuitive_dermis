function [tex_rs] = fill_array(tex_rs_t,tex_cnt,trialcnt,tdelaycnt,mr_cnt,scandir) 
% DEFAULT DIRECTION scandir=0 - skin sliding upward over stationary texture
% 1 2 3 ...
% 9 10 11 ...
mr_cnt_oneside = sqrt(mr_cnt) %% assumes array is square

    for k = 1:tex_cnt
        for trial=1:trialcnt
           disp(['Reading stimulus # ',num2str(k),' trial # ',num2str(trial)]);
           td = tdelaycnt; 
           if mr_cnt > 1
               for h=0:mr_cnt/mr_cnt_oneside - 1                
                    disp(['--->Filling receptors # ',num2str(mr_cnt_oneside*h+1),'-',num2str(mr_cnt_oneside*h+mr_cnt_oneside)]);
                    for i=1:mr_cnt_oneside
                        if tex_cnt==1
                            tex_rs(mr_cnt_oneside*h+i,:) = tex_rs_t{td-i+1};
                            if scandir==1
                                tex_rs(mr_cnt_oneside*h+i,:) = flip(tex_rs(mr_cnt_oneside*h+i,:));
                            end
                        else
                            tex_rs{k,trial}(mr_cnt_oneside*h+i,:) = tex_rs_t{k,trial,td-i+1};
                            if scandir==1
                                tex_rs{k,trial}(mr_cnt_oneside*h+i,:) = flip(tex_rs{k,trial}(mr_cnt_oneside*h+i,:));
                            end
                        end    
                    end
                    td = td-1;
               end
           else
               disp(['--->Filling LONE receptors']);
               if tex_cnt==1
                   tex_rs(1,:) = tex_rs_t{td};
                    if scandir==1
                        tex_rs(1,:) = flip(tex_rs(1,:));
                    end
               else
                   tex_rs{k,trial}(1,:) = tex_rs_t{k,trial,td};
                    if scandir==1
                        tex_rs{k,trial}(1,:) = flip(tex_rs{k,trial}(1,:));
                    end             
               end
           end
        end
    end

    size(tex_rs)
    if tex_cnt~=1
        size(tex_rs{1})
    end
end