function [p] = gx2_imhof_mod(x, w, k, lambda, s, m, RelTol, AbsTol, VPAflag)

    %% Modified version of "gx2_imhof.m"
    if size(x, 1) > 1
        error('size');
    end

    if VPAflag == true
        syms u
        imhof_integral = arrayfun(@(i) vpaintegral(@(u) gx2_imhof_integrand(u.', x(i), w', k', lambda(i), s, m, 'pdf'), ...
            u, 0, inf, 'AbsTol', AbsTol, 'RelTol', RelTol), 1:numel(x));
        p = imhof_integral / (2 * pi);
        if isa(p, 'sym')
            p = double(vpa(p));
        end
    else
        f_obj = @(u) gx2_imhof_integrand(u.', x, w', k', lambda, s, m, 'pdf');
        p = 1 / (2 * pi) * integral(f_obj, 0, inf, ArrayValued = true, AbsTol = AbsTol, RelTol = RelTol); %,
    end

end
