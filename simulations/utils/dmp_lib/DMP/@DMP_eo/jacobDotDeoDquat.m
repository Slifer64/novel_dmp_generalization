function dJ_deo_dQ = jacobDotDeoDquat(Qe, rotVel)

    deo = DMP_eo.rotVel2deo(rotVel, Qe);

    if ( (1-abs(Qe(1))) <= DMP_eo.zero_tol)
        dJ_deo_dQ = [-deo/3 zeros(3,3)];
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
    temp = (th*c_th-s_th)/s_th^2;

    dJ_deo_dQ = zeros(3,4);
    dJ_deo_dQ(:,1) = ((-th/s_th - 2*c_th*temp/s_th)*Eta + temp*(eye(3,3)-Eta)/th)*deo;
    dJ_deo_dQ(:,2:4) = -temp*dot(eta,deo)*eye(3,3);

%             dQe = -0.5*quatProd(Qe, [0; rotVel]);
%             dJ = getQToLogJacobianAcceleration( Qe, dQe );
%             dJ_deo_dQ2 = 2*dJ;
%             dJ_deo_dQ-dJ_deo_dQ2;
end
