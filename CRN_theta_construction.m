function theta = CRN_theta_construction(start_point, coord_cell, elementary_basis, F)
number_of_elementary_path = rank(elementary_basis);
stopover = sym(start_point);
theta = sym(1);
syms jj

for l = 1:number_of_elementary_path
    tmp = stopover + (jj-1) * elementary_basis(l,:)';
    stopover_cell = num2cell(tmp);
%    assumeAlso(jj <= coord_struct{l});
    theta = theta * symprod(F{l}(stopover_cell{:}), jj, 1, coord_cell{l});
    if l > 1
        stopover = stopover + coord_cell{l-1} * elementary_basis(l-1,:)';
    end
end
end

