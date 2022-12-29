function dJ_dQ_deo = jacobDotDquatDeo(Qe, rotVel)
            
    deo = DMP_eo.rotVel2deo(rotVel, Qe);
    
    if ( (1-abs(Qe(1))) <= DMP_eo.zero_tol)
        dJ_dQ_deo = [-deo'/4; zeros(3,3)];
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
    I_eta = eye(3,3) - Eta;
    temp = ((th*c_th-s_th)/th^2);
    
    dJ_dQ_deo = zeros(4,3);
    dJ_dQ_deo(1,:) = -0.25 * deo' * (c_th*Eta + (s_th/th)*I_eta);
    dJ_dQ_deo(2:4,:) = 0.25*dot(eta,deo)*( temp*I_eta - s_th*Eta ) + 0.25*temp*( eta*(deo'*I_eta) + (I_eta*deo)*eta' );

end
