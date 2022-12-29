function rotAccel = ddeo2rotAccel(ddeo, rotVel, Qe)
            
    deo = DMP_eo.rotVel2deo(rotVel, Qe);
    invQe = quatInv(Qe);
    J = DMP_eo.jacobDquatDeo(Qe);
    dJ = DMP_eo.jacobDotDquatDeo(Qe, rotVel);
    
    % rotVelQ = [0; rotVel];

    % rotAccel = -quatProd(quatProd(rotVelQ, invQe), J*deo) - 2*quatProd(invQe, dJ*deo) - 2*quatProd(invQe, J*ddeo);
    % rotAccel2 = 0.5*quatProd(rotVelQ, rotVelQ) - 2*quatProd(invQe, dJ_dQ_deo*deo) - 2*quatProd(invQe, J_dQ_deo*ddeo);
    rotAccel = - 2 * (quatProd(invQe, dJ*deo) + quatProd(invQe, J*ddeo));
    rotAccel = rotAccel(2:4);
    
end
