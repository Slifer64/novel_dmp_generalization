function deo = rotVel2deo(rotVel, Qe)
            
    J_deo_dQ = DMP_eo.jacobDeoDquat(Qe);
    deo = -0.5*J_deo_dQ * quatProd(Qe, [0; rotVel]);

end
