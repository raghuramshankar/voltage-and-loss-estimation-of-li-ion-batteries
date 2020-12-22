function objec = opt(x, V_batt_time, I, z, soc, dt, s, i_r, v_c, ind, v_t, ocv)
f = exp(-dt/(x(2)*x(3)));
    for k = 1:s-1
        i_r(k+1) = diag(f)*i_r(k) + (1-f)*I(k);
        v_c(k+1) = x(2).*i_r(k+1);
        ind(k) = interp1(soc, 1:length(soc), z(k), 'nearest');
        v_t(k+1) = ocv(ind(k)) - v_c(k) - I(k).*x(1);
    end
    objec = sqrt(mean((V_batt_time - v_t).^2));
end