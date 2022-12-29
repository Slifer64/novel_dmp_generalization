%% N-DoF GMP Update class
%

classdef GMPo_Update < GMP_Update

    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMPo_Update(gmp)
                
            if ( ~ strcmpi(class(gmp),'GMPo') )
                error(['[GMPo_Update::GMPo_Update]: argument of type ''GMPo'' was expected, but a type ''' class(gmp) ''' was provided instead']);
            end
                
            this@GMP_Update(gmp);

        end

        % function initSigmaw(this)
        
        % function initSigmaWfromMsr(this, x_data)
        
        % function enableSigmawUpdate(this, flag)
        
        % function setSigmaW(this, Sw)
        
        % function Sw = getSigmaW(this)
        
        % function setMsrNoiseVar(this, rv)
        
        
        %% ==================================================
        %% ===============   Online update  =================
        %% ==================================================
        
        % function updatePos(this, x, y, r_n)
        
        % function updateVel(this, x, x_dot, y_dot, r_n)
        
        % function updateAccel(this, x, x_dot, x_ddot, y_ddot, r_n)
        
        function updateQuat(this, x, Q, r_n)
            
            if (nargin < 4), r_n=this.rv; end
            
            y = GMPo.quat2q(Q, this.gmp.getQ0());
            this.updatePos(x, y, r_n);
            
        end
        
        function updateRotVel(this, x, x_dot, rotVel, Q, r_n)
            
            if (nargin < 6), r_n=this.rv; end
            
            Q1 = GMPo.getQ1(Q, this.gmp.getQ0());
            y_dot = gmp_.rotVel_to_qLogDot(rotVel, Q1);
            this.updateVel(x, x_dot, y_dot, r_n);
            
        end
        
        function updateRotAccel(this, x, x_dot, x_ddot, rotAccel, rotVel, Q, r_n)
            
            if (nargin < 8), r_n=this.rv; end
            
            Q1 = GMPo.getQ1(Q, this.gmp.getQ0());
            y_ddot = gmp_.rotAccel_to_qLogDDot(rotAccel, rotVel, Q1);
            this.updateAccel(x,x_dot,x_ddot, y_ddot, r_n);
          
        end

        % function updateWeights(this, s, Z, type, r_n)
        
        % function plotWeightsCovariance(this, ax)
    end
    
end
