function objec = opt(x, r0, V_batt_time, I, z, soc, dt, s, v_t, i_r, v_c, ocv, n)
    for k = 1:n
        r(k) = x(2*k-1);
        c(k) = x(2*k);
    end
    for j = 1:n
        f(j) = exp(-dt/(r(j)*c(j)));
    end
    for k = 1:s-1
        i_r(:, k+1) = diag(f)*i_r(:, k) + (1-f')*I(k);
        v_c(:, k) = i_r(:, k).*r';
        ind = interp1(soc, 1:length(soc), z(k), 'nearest');
        v_t(k+1) = ocv(ind) - sum(v_c(:, k)) - I(k).*r0;
    end
    v_t = reshape(v_t, size(V_batt_time));
    objec = sqrt(mean((V_batt_time - v_t).^2));
end