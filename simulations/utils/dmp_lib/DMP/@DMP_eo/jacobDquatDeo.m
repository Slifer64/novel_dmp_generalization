function J_dQ_deo = jacobDquatDeo(Qe)
            
    if ( (1-abs(Qe(1))) <= DMP_eo.zero_tol)
        J_dQ_deo = [zeros(1, 3); eye(3,3)];
        return;
    end

    w = Qe(1);
    v = Qe(2:4);
    norm_v = norm(v);
    eta = v / norm_v;
    s_th = norm_v;
    c_th = w;
    th = atan2(s_th, c_th);
    Eta = eta*eta';
    
    J_dQ_deo = zeros(4,3);
    J_dQ_deo(1,:) = -0.5 * s_th * eta';
    J_dQ_deo(2:4,:) = 0.5 * ( (eye(3,3) - Eta)*s_th/th + c_th*Eta );

end
