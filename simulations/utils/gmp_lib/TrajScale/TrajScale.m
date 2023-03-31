%% Trajectory scaling class
%  Calculates the scaling matrix for scaling spatially a trajectory from a  
%  nominal start/final position {Y0d, Ygd} to a new start/final position {Y0, Yg}.
%  The scaling is performed according to the scaling method (see @TrajScale::ScaleMethod enum)
%  The scaling is returned as matrix T_sc that premultiplies a position to
%  scale it, i.e. Y_scaled = T_sc * Y
%  The scaling is also invertible, i.e. Y = inv(T_sc) * Y_scaled

classdef TrajScale < matlab.mixin.Copyable
    
    properties (Constant, Access = public)
        
        % enum TrajScaleType{
        NONE = -1
        PROP_SCALE = 0
        ROT_MIN_SCALE = 1
        ROT_WB_SCALE = 2
        ROT_DIAL_SCALE = 3
        %}
        
    end
    
    methods (Access = public)
        
        %% Constructor.
        %  @param[in] n_dof: degrees of freedom.
        function this = TrajScale(n_dof)
            
            this.n_dof = n_dof;
            this.Y0d = zeros(n_dof,1);
            this.Ygd = ones(n_dof,1);
            % Initialization required for 'setScaleMethod'
            % will be overwirtten then on the 'setInitFinalPos' call
            this.Y0 = zeros(n_dof,1); 
            this.Yg = ones(n_dof,1);
            
            this.setNominalStartFinalPos(this.Y0d, this.Ygd);
            this.setNewStartFinalPos(this.Y0, this.Yg);
            
        end
        
        %% Returns the matrix that scales to the new init/target position, according to the set scaling method.
        %  @return: The scaling matrix to the new init/target position.
        function T_sc = getScaling(this)
            
            T_sc = this.T_sc;
            
        end
        
        %% Returns the inverse of the matrix that scales to the new init/target position, according to the set scaling method.
        %  @return: The inverse scaling matrix to the new init/target position.
        function inv_T_sc = getInvScaling(this)
            
            inv_T_sc = this.inv_T_sc;
            
        end

        %% Set the new start and final position.
        %  @param[in] Y0: new start position.
        %  @param[in] Yg: new final position.
        function setNewStartFinalPos(this, Y0, Yg)

            this.T_sc = this.calcScaling(Y0, Yg);
            this.inv_T_sc = this.calcInvScaling();
            
        end
        
        %% Set the nominal start and final position.
        %  @param[in] Y0d: nominal start position.
        %  @param[in] Ygd: nominal final position.
        function setNominalStartFinalPos(this, Y0d, Ygd)
            
            this.Y0d = Y0d;
            this.Ygd = Ygd;
            
            this.T_sc = this.calcScaling(this.Y0, this.Yg);
            this.inv_T_sc = this.calcInvScaling();
            
        end
        
        function n_dofs = numOfDoFs(this)
            
            n_dofs = length(this.Y0d);
            
        end
        
        function cp_obj = deepCopy(this)
           
            cp_obj = this.copy();
            
        end
        
        function y0 = getY0(this) 
            y0 = this.Y0;
        end
        
        function yg = getYg(this) 
            yg = this.Yg;
        end
        
    end
    
    methods (Abstract, Access = public)
        
        scale_type = getScaleType(this)
        
    end
    
    methods (Abstract, Access = protected)
        
        sc = calcScaling(this)
        
        inv_sc = calcInvScaling(this)
        
    end
    
    properties (Access = protected)
        
        n_dof % number of DoFs
        
        Y0d % Nominal start position.
        Ygd % Nominal final position.
        Y0 % New start position.
        Yg % New final position.
        
        T_sc % scaling matrix, calculated based on the set scale method
        inv_T_sc % inverse scaling matrix, calculated based on the set scale method

    end
    
end
