function [pop,pop_NA] =  gene_layout_by_indices_one(wf,indices)

    pop_NA = zeros(1, wf.rows*  wf.cols);
    pop = pop_NA;
    pop(indices) = 1;
    pop_NA(indices) = 1;
    pop_NA(wf.NA_loc)= 2;

end