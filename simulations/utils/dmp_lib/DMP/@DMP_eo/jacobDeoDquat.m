function J_deo_dQ = jacobDeoDquat(Qe)
            
    if ( (1-abs(Qe(1))) <= DMP_eo.zero_tol)
        J_deo_dQ = [zeros(3,1) eye(3,3)];
        return;
    end
    
    w = Qe(1);
    v = Qe(2:4);
    norm_v = norm(v);
    eta = v / norm_v;
    s_th = norm_v;
    c_th = w;
    th = atan2(s_th, c_th);
    
    J_deo_dQ = zeros(3,4);
    J_deo_dQ(:,1) = 2*eta*(th*c_th - s_th)/s_th^2;
    J_deo_dQ(:,2:4) = 2*eye(3,3)*th/s_th;

end
