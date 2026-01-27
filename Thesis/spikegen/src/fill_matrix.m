function [tex_out] = fill_matrix(tex_in,inds,t)

    ti = ismember(t,inds);
    for t = 1:length(ti)
        if ti(t) == 1
            tex_out(t) = tex_in(t);
        else
            tex_out(t) = 0;
        end
    end

end