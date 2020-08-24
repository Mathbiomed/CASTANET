function theta = CRN_theta_construction(start_point, coord_struct, elementary_basis, F)
number_of_elementary_path = rank(elementary_basis);
stopover = sym(start_point);
theta = sym(1);
syms jj

for l = 1:number_of_elementary_path
    tmp = stopover + (jj-1) * elementary_basis(l,:);
    stopover_cell = num2cell(tmp);
    theta = theta * symprod(F{l}(stopover_cell{:}), jj, 1, coord_struct{l});
    if l > 1
        stopover = stopover + coord_struct{l-1} * elementary_basis(l-1,:);
    end
end
end

