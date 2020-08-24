function path_final = CRN_find_path(b,complexes)
% complexes: d * |C| matrix whose columns are complex vectors. 
% b : reference (or start) point vector.

[d, num_complex] = size(complexes);
num_unit_path = num_complex * (num_complex - 1) / 2;
unit_path = zeros(d, num_unit_path);
ai = 1; bi = 1;
coordinate = zeros(2, num_unit_path);
for ii = 1:num_unit_path
    ai = ai + 1;
    if ai > num_complex
        bi = bi + 1;
        ai = bi + 1;
    end
    unit_path(:,ii) = complexes(:, ai) - complexes(:, bi);
    coordinate(:, ii) = [ai; bi];
end

% This symbolic computation part is needed to be edited. 
% This part might not be able to be properly coded.
n = sym('n', [d 1], 'integer');
X = sym('x', [num_unit_path, 1], 'integer');
eqn = unit_path * X == n - b;
path_index = solve(eqn, X);

path_final= [];
for ii = 1:num_unit_path
    if path_index(ii) > 0
        for jj = 1:path_index(ii)
            path_final = [path_final, coordinate(:,ii)];
        end
    elseif path_index(ii) < 0
        for jj = 1:(-path_index(ii))
            path_final = [path_final, flip(coordinate(:,ii))];
        end
    end
end
% path_final: 2 * m matrix contaning a chain from b to n of length m.
% each column is [aj, bj] which means (nu_aj - nu_bj) is a unit path
% consisting the whole path.

% NEED editing to perform symbolic computation properly. %%%%%%%%
end

