function [py, preC, logE] = spa_pdf_vec(y, nuW, p1, d, nuN2)

    %SPA: https://infoscience.epfl.ch/server/api/core/bitstreams/3154795e-11b3-4e0d-98c7-70cac92c2144/content
    % Saddle-point approximation (SPA) for generalized non-central chi-square (gXi2)
    % - observations: y
    % -
    % - gXi2 order: d

    y = y(:); %COLUMN
    p1 = p1(:); %COLUMN

    logMGF = @(t, w, d, lambda) ...
        (t .* w .* lambda) ./ (1 - 2 * w * t) + t.^2 * nuN2 / 2 - d / 2 * log(1 - 2 * w * t);

    logMGFd1 = @(t, w, d, lambda) ...
        (w * (lambda + d)) ./ (1 - 2 * w * t) + t * nuN2 + ...
        + (2 * w.^2 * lambda * t) ./ ((1 - 2 * w * t).^2);

    n = length(y);

    %Don't do this by hand... use MATLAB symbolic computation
    %syms y w lambda t nuN2 z d
    %Define: eqn = logMGFd1 - y == 0
    %coeffs(solve(eqn))
    logMGFd1_Y_poly_coeff = @(y, w, d, lambda) [4 * nuN2 * w^2 * ones(n, 1), - 4 * y * w^2 - 4 * nuN2 * w, - 2 * d * w^2 + 4 * y * w + nuN2, d * w - y + lambda * w];

    logMGFd2 = @(t, w, d, lambda) + (2 * w^2 * (2 * lambda + d)) ./ (1 - 2 * w * t).^2 + nuN2 ...
        + 8 * (w.^3 * lambda .* t) ./ ((1 - 2 * w * t).^3);

    w = nuW / 2;
    lamb = abs(p1).^2 / w;

    % Methods:
    % SPA Cardano Vectorized
    mypolMAT = logMGFd1_Y_poly_coeff(y, w, d, lamb); % n x (polyorder + 1)
    t_hat = wdenoise.cardano_roots_vec(mypolMAT, 1 / (2 * w) - eps);


    % Outputs:
    preC = (1 ./ sqrt(2 * pi * logMGFd2(t_hat, w, d, lamb))).';
    logE = (logMGF(t_hat, w, d, lamb) - t_hat .* y).';
    py = (preC .* exp(logE)).';

end
