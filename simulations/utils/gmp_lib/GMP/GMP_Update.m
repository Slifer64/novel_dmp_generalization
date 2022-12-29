%% N-DoF GMP Update class
%

classdef GMP_Update < matlab.mixin.Copyable

    methods (Access = public)
        
        %% GMP constructor.
        %  @param[in] gmp: n_DoF dmp.
        function this = GMP_Update(gmp)
                
%             if ( ~ strcmpi(class(gmp),'GMP') )
%                 error(['[GMP_Update::GMP_Update]: argument of type ''GMP'' was expected, but a type ''' class(gmp) ''' was provided instead']);
%             end
            
            this.gmp = gmp;
            
            this.initSigmaw();
            this.rv = 1;
            
            this.recursiveUpdate(false);
            this.syncUpdate(true);

        end

        function initSigmaw(this)
            
            N_kernels = this.gmp.numOfKernels();
            this.Sigma_w = eye(N_kernels, N_kernels);
            
        end
        
        function initExpSigmaw(this, decay_rate)
            
            if (nargin < 2), decay_rate = 0.01; end
            
            N_kernels = this.gmp.numOfKernels();
            S = zeros(N_kernels,N_kernels);
            for i=1:N_kernels
                for j=i+1:N_kernels
                    S(i,j) = exp(-decay_rate * abs(i-j));
                end
            end
            this.Sigma_w = S + S' + eye(N_kernels,N_kernels); 
  
        end
        
        function initSigmaWfromMsr(this, x_data, tol)
              
            if (nargin < 3), tol=1e-2; end
            
            n_data = length(x_data);
            H = zeros(this.gmp.numOfKernels(), n_data);
            for j=1:n_data, H(:,j) = this.gmp.regressVec(x_data(j)); end
            
            this.Sigma_w = inv(H*H' + tol*eye(this.gmp.numOfKernels()));
            
        end
        
        function initSigmaWfromVelMsr(this, x_data, tau, tol)
            
            if (nargin < 4), tol=1e-2; end
            
            n_data = length(x_data);
            H = zeros(this.gmp.numOfKernels(), n_data);
            for j=1:n_data, H(:,j) = this.gmp.regressVecDot(x_data(j), 1/tau); end
            
            this.Sigma_w = inv(H*H' + tol*eye(this.gmp.numOfKernels()));
            
        end
        
        function initSigmaWfromAccelMsr(this, x_data, tau, tol)
            
            if (nargin < 4), tol=1e-2; end
            
            n_data = length(x_data);
            H = zeros(this.gmp.numOfKernels(), n_data);
            for j=1:n_data, H(:,j) = this.gmp.regressVecDDot(x_data(j), 1/tau, 0); end
            
            this.Sigma_w = inv(H*H' + tol*eye(this.gmp.numOfKernels()));
            
        end
        
        function recursiveUpdate(this, set)
           
            this.recursive_up = set;
            
        end
        
        function syncUpdate(this, set)
           
            this.sync_up = set;
            
        end
        
        function setSigmaW(this, Sw)
        
            this.Sigma_w = Sw;
            
        end
        
        function Sw = getSigmaW(this)
           
            Sw = this.Sigma_w;
            
        end
        
        function setMsrNoiseVar(this, rv)
           
            this.rv = rv;
            
        end
        
        %% ==================================================
        %% ===============   Online update  =================
        %% ==================================================
        
        function updateNow(this)
            
            this.updateWeights(this.batch_s, this.batch_Z, this.batch_type, this.batch_r_n);
            this.clearBatch();
    
        end

        function updatePos(this, x, y, r_n)
            
            if (nargin < 4), r_n=this.rv; end
            
            this.batch_s = [this.batch_s, GMP_phase(x,0,0)];
            this.batch_Z = [this.batch_Z, y];
            this.batch_type = [this.batch_type, GMP_UpdateType.POS];
            this.batch_r_n = [this.batch_r_n, r_n];

            if (this.sync_up), this.updateNow(); end
            
        end
        
        function updateVel(this, x, x_dot, y_dot, r_n)
            
            if (nargin < 5), r_n=this.rv; end

            this.batch_s = [this.batch_s, GMP_phase(x,x_dot,0)];
            this.batch_Z = [this.batch_Z, y_dot];
            this.batch_type = [this.batch_type, GMP_UpdateType.VEL];
            this.batch_r_n = [this.batch_r_n, r_n];

            if (this.sync_up), this.updateNow(); end
            
        end
        
        function updateAccel(this, x, x_dot, x_ddot, y_ddot, r_n)
            
            if (nargin < 6), r_n=this.rv; end

            this.batch_s = [this.batch_s, GMP_phase(x,x_dot,x_ddot)];
            this.batch_Z = [this.batch_Z, y_ddot];
            this.batch_type = [this.batch_type, GMP_UpdateType.ACCEL];
            this.batch_r_n = [this.batch_r_n, r_n];

            if (this.sync_up), this.updateNow(); end
            
        end
        
        
        %% Updates the weights based on 'n' measurements.
        %  @param[in] s: 1 x n vector of type GMP_phase, where the j-th entry is the phase for the j-th measurement.
        %  @param[in] Z: n_dof x n matrix, where the j-th column is the j-th measurement.
        %  @param[in] type: 1 x n vector of type GMP_UpdateType, where the j-th entry has the update type for the j-th measurement.
        %  @param[in] r_n: 1 x n vector, where the j-th entry is the noise variance for the j-th measurement (optional, default=this.rv).
        function updateWeights(this, s, Z, type, r_n)
 
            if (nargin < 5), r_n = this.rv; end
                
            n = size(Z,2);
            n_ker = this.gmp.numOfKernels();
            sc = this.gmp.getScaling();
            inv_sc = this.gmp.getInvScaling();
            
            if (length(r_n) == 1), r_n = repmat(r_n(1), 1,n); end
            Rn = diag(r_n);
 
            H = zeros(n_ker, n);
            
            for j=1:n
                if (type(j) == GMP_UpdateType.POS)
                    H(:,j) = this.gmp.regressVec(s(j).x);
                    Z(:,j) = Z(:,j) - this.gmp.Y0 + sc*this.gmp.Y0d;
                elseif (type(j) == GMP_UpdateType.VEL)
                    H(:,j) = this.gmp.regressVecDot(s(j).x, s(j).x_dot);
                else % (type == GMP_UpdateType.ACCEL)
                    H(:,j) = this.gmp.regressVecDDot(s(j).x, s(j).x_dot, s(j).x_ddot);
                end
                %Z(:,j) = Z(:,j)./sc;
            end
            Z = inv_sc*Z;
            
            C = this.Sigma_w*H;
            dS = decomposition(H'*this.Sigma_w*H + Rn, 'ldl', 'upper');
            B = dS \ C';
            % B = eye(n,n) / (H'*this.Sigma_w*H + Rn) * C';
            this.gmp.W = this.gmp.W + (Z - this.gmp.W*H)*B;
            
            if (this.recursive_up)
                % this.Sigma_w = this.Sigma_w - C*B;
                
                % use "Joseph form" for numerical stability
                I_HB = eye(size(this.Sigma_w)) - H*B;
                this.Sigma_w = I_HB'*this.Sigma_w*I_HB + B'*Rn*B;

                % to ensure symmetry
                this.Sigma_w  = (this.Sigma_w + this.Sigma_w') / 2;
            end
            
        end

        function plotWeightsCovariance(this, ax)
            
            if (nargin < 2)
                figure;
                ax = axes();
            end
            imshow(this.Sigma_w, [min(this.Sigma_w(:)) max(this.Sigma_w(:))], 'InitialMagnification', 2000, 'Parent',ax);
            
        end
        
    end
    
    methods (Access = protected)
       
        function clearBatch(this)
            
            this.batch_s = [];
            this.batch_Z = [];
            this.batch_type = [];
            this.batch_r_n = [];
            
        end
        
    end
    
    properties (Access = protected)

        gmp % n_DoF GMP
        
        Sigma_w % weights covariance for trajectory update
        rv % noise variance for trajectory update
        recursive_up
        sync_up
        
        batch_s
        batch_Z
        batch_type
        batch_r_n

    end
    
end
