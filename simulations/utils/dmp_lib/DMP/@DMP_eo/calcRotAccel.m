function rotAccel = calcRotAccel(this, x, Q, rotVel, Qg, tau_dot, yc, zc, yc_dot)
            
    if (nargin < 6), tau_dot = 0; end
    if (nargin < 7), yc = zeros(3,1); end
    if (nargin < 8), zc = zeros(3,1); end
    if (nargin < 9), yc_dot = zeros(3,1); end

    a_z = [this.dmp{1}.a_z; this.dmp{2}.a_z; this.dmp{3}.a_z];
    b_z = [this.dmp{1}.b_z; this.dmp{2}.b_z; this.dmp{3}.b_z];
    tau = this.getTau();
    
    
    Qg_prev = this.Qg;
    this.setQg(Qg);

    Qe = DMP_eo.quatError(Q, Qg);
    eo = DMP_eo.quat2eo(Q, Qg);
    invQe = quatInv(Qe);

    rotVelQ = [0; rotVel];
    QeRotVel = quatProd(Qe,rotVelQ);

    J_deo_dQ = DMP_eo.jacobDeoDquat(Qe);
    J_dQ_deo = DMP_eo.jacobDquatDeo(Qe);
    dJ_dQ_deo = DMP_eo.jacobDotDquatDeo(Qe, rotVel);

    fo = zeros(3,1);
    for i=1:3, fo(i) = this.dmp{i}.shapeAttractor(x, 0); end

    deo = DMP_eo.rotVel2deo(rotVel, Qe);
    ddeo = (-a_z.*b_z.*eo - tau*a_z.*deo + a_z.*yc + fo + tau*yc_dot - tau*tau_dot*deo + zc)/tau^2;

    rotAccel1 = quatProd(invQe, dJ_dQ_deo*J_deo_dQ*QeRotVel);
    % rotAccel2 = 2*quatProd(invQe, J_dQ_deo* (a_z.*b_z.*eo - 0.5*tau*a_z.*(J_deo_dQ*QeRotVel) - fo) ) / tau^2;
    rotAccel2 = 2*quatProd(invQe, J_dQ_deo*-ddeo);   
    rotAccel = rotAccel1 + rotAccel2;  
    rotAccel = rotAccel(2:4);
    
    this.setQg(Qg_prev);

	% norm(rotAccel-rotAccel_2)

end