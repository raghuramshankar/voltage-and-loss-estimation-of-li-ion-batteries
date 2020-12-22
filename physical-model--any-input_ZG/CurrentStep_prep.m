
% initialize current distribution
for i = 1:5
    R_x_neg = spdiags([-R_ct_neg(2:end)' (R_s_neg + R_ct_neg(1:end-1)+R_ct_neg(2:end)+R_l_neg)' -R_ct_neg(1:end-1)'],-1:1,n_mesh_neg,n_mesh_neg);
    V_x_neg = phis_U_neg(1:end-1)-phis_U_neg(2:end)+I*R_l_neg;
    V_x_neg(1) = V_x_neg(1)+I*R_ct_neg(1);
    I_x_neg = R_x_neg\V_x_neg';
    
    I_ct_neg = [-diff([I;I_x_neg]);I_x_neg(end)]'; % charge transfer current per electrode area
    I_loc_neg = I_ct_neg./(Sa_neg*volume_length_neg); % charge transfer current per active area
    i_0_update = F*k_ct_neg*(c_s_max_neg-c_s_surface_neg).^0.5.*c_s_surface_neg.^0.5.*c_l(1:n_mesh_neg+1).^0.5;
    eta_neg = 2*R*T/F*asinh(I_loc_neg./2./i_0_update);
    R_ct_neg = eta_neg./I_ct_neg;
end
%==========================================================================


% initialize current distribution
for i = 1:5
    R_x_pos = spdiags([-R_ct_pos(2:end)' (R_s_pos + R_ct_pos(1:end-1)+R_ct_pos(2:end)+R_l_pos)' -R_ct_pos(1:end-1)'],-1:1,n_mesh_pos,n_mesh_pos);
    V_x_pos = phis_U_pos(1:end-1)-phis_U_pos(2:end)+I*R_l_pos;
    V_x_pos(end) = V_x_pos(end)+I*R_ct_pos(end);
    I_x_pos = R_x_pos\V_x_pos';
    
    I_ct_pos = [I_x_pos(1);diff([I_x_pos;I])]'; % charge transfer current per electrode area
    I_loc_pos = I_ct_pos./(Sa_pos*volume_length_pos); % charge transfer current per active area
    i_0_update = F*k_ct_pos*(c_s_max_pos-c_s_surface_pos).^0.5.*c_s_surface_pos.^0.5.*c_l(n_mesh_neg+n_mesh_sep+3:end).^0.5;
    eta_pos = 2*R*T/F*asinh(I_loc_pos./2./i_0_update);
    R_ct_pos = eta_pos./I_ct_pos;
end


