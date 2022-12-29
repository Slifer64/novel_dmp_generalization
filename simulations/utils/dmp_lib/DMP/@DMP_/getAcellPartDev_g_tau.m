function dC_dtheta = getAcellPartDev_g_tau(this, t, y, dy, y0, x, g, tau)

    dC_dtheta = zeros(2,1);

    K_dmp = this.a_z*this.b_z;
    D_dmp = this.a_z;
    psi = this.kernelFunction(x);
    sum_psi = sum(psi) + this.zero_tol;
    sum_w_psi = psi'*this.w;
    shape_attr_gat = this.shapeAttrGating(x);

    theta1 = g;
    theta2 = 1/tau;

    dshape_attr_gat_dtheta2 = this.shape_attr_gating_ptr.getPartDev_1oTau(t,x);

    dPsidtheta2 = -2*t*this.h.*(theta2*t-this.c).*psi;
    sum_w_dPsidtheta2 = this.w'*dPsidtheta2;
    dSumWPsi_dtheta2 = (sum_w_dPsidtheta2*sum_psi - sum_w_psi*sum(dPsidtheta2) ) / sum_psi^2;

    dC_dtheta(1) = (K_dmp + shape_attr_gat*sum_w_psi/sum_psi)*theta2^2;

    dC_dtheta(2) = 2*theta2* (K_dmp*(theta1-y) + shape_attr_gat*(theta1-y0)*sum_w_psi/sum_psi) + ...
        -D_dmp*dy + theta2^2*(theta1-y0)*( dshape_attr_gat_dtheta2*sum_w_psi/sum_psi + shape_attr_gat*dSumWPsi_dtheta2 );
    dC_dtheta(2) = dC_dtheta(2)*(-1/tau^2);

end