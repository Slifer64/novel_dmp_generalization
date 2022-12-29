function ddeo = rotAccel2ddeo(rotAccel, rotVel, Qe)
            
    rotVelQ = [0; rotVel];
    rotAccelQ = [0; rotAccel];
    
    J = DMP_eo.jacobDeoDquat(Qe);
    dJ = DMP_eo.jacobDotDeoDquat(Qe, rotVel);

    ddeo = -0.5 * (dJ * quatProd(Qe, rotVelQ) + J * quatProd( Qe, rotAccelQ-0.5*quatProd(rotVelQ,rotVelQ) ) );

end
