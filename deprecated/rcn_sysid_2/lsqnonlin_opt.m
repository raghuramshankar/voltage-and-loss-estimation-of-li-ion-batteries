function objec = lsqnonlin_opt(x, V_batt_time_rc, I_rc, dt, s, i_r, n, i)
    for k = 1:n
        r(k) = x(2*k-1);
        c(k) = x(2*k);
    end
    for j = 1:n
        f(j) = exp(-dt/(r(j).*c(j)));
    end
    for k = 1:s-i-3
        i_r(:, k+1) = diag(f)*i_r(:, k) + (ones(n, 1)-f')*I_rc(k);
        v_c(:, k+1) = i_r(:,k).*r';
        v_t(k) = sum(v_c(:, k));
    end
%     objec = sqrt(mean((V_batt_time_rc - v_t').^2));
    objec = V_batt_time_rc - v_t;
end