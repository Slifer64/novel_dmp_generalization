function rotVel = deo2rotVel(deo, Qe)
            
    J_dQ_deo = DMP_eo.jacobDquatDeo(Qe);
    rotVel = -2 * quatProd( quatInv(Qe), J_dQ_deo*deo );
    rotVel = rotVel(2:4);
    
end
