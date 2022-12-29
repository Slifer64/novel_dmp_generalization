%% DMP class
%  Implements an 1-D this.
%  The DMP's evolution is driven the phase varialbe 'x',  which is produced
%  by a linear canonical clock. The phase variable goes from x(0)=0 to x(tau)=1,
%  where 'tau' is the total movement's duration.
%  The DMP has the in general the following form:
%
%     tau*dz = g1(x)*( a_z*(b_z*(g-y) - z ) + g2(x)*fs*f(x) + z_c
%     tau*dy = z + y_c;
%
%  where
%     tau: is scaling factor defining the duration of the motion
%     a_z, b_z: constants relating to a spring-damper system
%     fs: scaling of the forcing term (typically fs = g0-y0)
%     g: the goal-final position
%     y0: the initial position
%     x: the phase variable
%     y,dy,ddy: the position, velocity and accelaration of the motion
%     f(x): the forcing term defined by the normalized weighted sum of the 
%        kernel functions (gaussian kernels), i.e.:
%        f(x) = w'*Psi(x)/ sum(Psi(x));
%     g1(x): the gating factor of the linear term
%     g2(x): the gating factor of the non-linear forcing term
%

classdef DMP_ < matlab.mixin.Copyable

    methods (Access = public)
        
        %% DMP constructor.
        %  @param[in] N_kernels: the number of kernels
        %  @param[in] a_z: Parameter 'a_z' relating to the spring-damper system.
        %  @param[in] b_z: Parameter 'b_z' relating to the spring-damper system.
        %  @param[in] can_clock_ptr: Pointer to a DMP canonical system object.
        function this = DMP_(N_kernels, a_z, b_z, can_clock_ptr, shape_attr_gating_ptr)

            % if (nargin < 4), can_clock_ptr = CanonicalClock(); end
            % if (nargin < 5), shape_attr_gating_ptr=SigmoidGatingFunction(1.0, 0.5); end
                
            this.N_kernels = N_kernels;
            this.a_z = a_z;
            this.b_z = b_z;
            this.can_clock_ptr = can_clock_ptr;
            this.shape_attr_gating_ptr = shape_attr_gating_ptr;
            
            this.zero_tol = 1e-30; %realmin;

            this.w = zeros(this.N_kernels,1);
            this.setCenters();
            this.setStds(1.0);
            this.setY0(0.0);
            
            this.dy = 0;
            this.dz = 0;
            this.dx = 0;
            
        end

        
        %% Trains the DMP.
        %  @param[in] Time: Row vector with the timestamps of the training data points.
        %  @param[in] yd_data: Row vector with the desired potition.
        %  @param[in] dyd_data: Row vector with the desired velocity.
        %  @param[in] ddyd_data: Row vector with the desired accelaration.
        %  @param[in] y0: Initial position.
        %  @param[in] g: Target-goal position.
        %
        %  \note The timestamps in \a Time and the corresponding position,
        %  velocity and acceleration data in \a yd_data, \a dyd_data and \a
        %  ddyd_data need not be sequantial in time.
        [train_error, F, Fd] = train(this, train_method, Time, yd_data, dyd_data, ddyd_data)
        
        
        %% Returns the derivatives of the DMP states.
        %  @param[in] x: phase variable.
        %  @param[in] y: \a y state of the this.
        %  @param[in] z: \a z state of the this.
        %  @param[in] y0: initial position.
        %  @param[in] g: Goal position.
        %  @param[in] y_c: coupling term for the dynamical equation of the \a y state.
        %  @param[in] z_c: coupling term for the dynamical equation of the \a z state.
        %  @param[out] dy: derivative of the \a y state of the this.
        %  @param[out] dz: derivative of the \a z state of the this.
        %  @param[out] dx: derivative of the phase variable of the this.
        update(this, x, y, z, g, y_c, z_c)
        
        
        function dx = getXdot(this), dx=this.dx; end
        function dy = getYdot(this), dy=this.dy; end
        function dz = getZdot(this), dz=this.dz; end
        
        
        %% Returns the DMP's acceleration.
        ddy = getYddot(this, tau_dot, yc_dot)
        
        
        ddy = calcYddot(this, x, y, dy, g, tau_dot, yc, zc, yc_dot)
        
        
        function n_ker = numOfKernels(this), n_ker = length(this.w); end

        
        function setY0(this, y0), this.y0 = y0; end

        
        %% Returns the time scale of the DMP.
        %  @param[out] tau: The time scale of the this.
        function tau = getTau(this), tau = this.can_clock_ptr.getTau(); end
        
        
        %% Sets the time scale of the DMP.
        function setTau(this, tau), this.can_clock_ptr.setTau(tau); end
   
        
        %% Returns the phase variable.
        %  @param[in] t: The time instant.
        %  @param[out] x: The phase variable for time 't'.
        function x = phase(this, t), x = this.can_clock_ptr.getPhase(t); end
   
        
        %% Returns the derivative of the phase variable.
        %  @param[in] x: The phase variable.
        %  @param[out] dx: The derivative of the phase variable.
        function dx = phaseDot(this, x), dx = this.can_clock_ptr.getPhaseDot(x); end

        
        %% Returns the goal attractor of the this.
        %  @param[in] x: The phase variable.
        %  @param[in] y: \a y state of the this.
        %  @param[in] z: \a z state of the this.
        %  @param[in] g: Goal position.
        %  @param[out] goal_attr: The goal attractor of the this.
        goal_attr = goalAttractor(this, x, y, z, g)

        
        %% Returns the partial derivative of the DMP's acceleration wrt to the goal and tau.
        %  @param[in] t: current timestamp.
        %  @param[in] y: position.
        %  @param[in] dy: velocity.
        %  @param[in] y0: initial position.
        %  @param[in] x_hat: phase variable estimate.
        %  @param[in] g_hat: goal estimate.
        %  @param[in] tau_hat: time scale estimate.
        %  @param[out] dC_dtheta: partial derivative of the DMP's acceleration wrt to the goal and tau.
        dC_dtheta = getAcellPartDev_g_tau(this, t, y, dy, y0, x, g, tau)

        
        %% Creates a deep copy of this object
        cp_obj = deepCopy(this)
 
    end
    
    methods (Access = protected)
        
        %% Returns the shape attractor gating factor.
        %  @param[in] x: The phase variable.
        function sAttrGat = shapeAttrGating(this, x)

            sAttrGat = this.shape_attr_gating_ptr.getOutput(x);
            sAttrGat(sAttrGat<0) = 0.0;

        end
        
        
        %% Returns the goal attractor gating factor.
        %  @param[in] x: The phase variable.
        function gAttrGat = goalAttrGating(this, x), gAttrGat = 1.0; end

        
        %% Returns the forcing term of the dmp.
        %  @param[in] x: The phase variable.
        %  @param[out] f: The normalized weighted sum of Gaussians.
        f = forcingTerm(this,x)
        
        
        %% Returns a column vector with the values of the kernel functions of the DMP.
        %  @param[in] x: phase variable
        %  @param[out] psi: column vector with the values of the kernel functions of the DMP
        psi = kernelFunction(this,x)
        
        
        %% Sets the centers for the kernel functions of the DMP according to the canonical system.
        setCenters(this)

        
        %% Sets the standard deviations for the kernel functions  of the DMP.
        %  Sets the variance of each kernel equal to squared difference between the current and the next kernel.
        %  @param[in] kernelStdScaling: Scales the std of each kernel by 'kernelStdScaling' (optional, default = 1.0).
        setStds(this, kernelStdScaling)
        
    end

    methods (Abstract, Access = public)
        
        shape_attr = shapeAttractor(this, x, g)

    end
    
    methods (Abstract, Access = protected)
        
        Fd = calcFd(this, x, y, dy, ddy, g)
        
        Fd = calcLearnedFd(this, x, g)
        
        f_scale = forcingTermScaling(this, g)

    end
    
    properties (Access = public)
        
        N_kernels % number of kernels (basis functions)

        a_z % parameter 'a_z' relating to the spring-damper system
        b_z % parameter 'b_z' relating to the spring-damper system

        can_clock_ptr % pointer to the canonical clock
        shape_attr_gating_ptr % pointer to gating function for the shape attractor
        
        w % N_kernelsx1 vector with the weights of the DMP
        c % N_kernelsx1 vector with the kernel centers of the DMP
        h % N_kernelsx1 vector with the kernel stds of the DMP

    end
    
    
    properties (Access = protected)

        zero_tol % tolerance value used to avoid divisions with very small numbers
        
        y0 % initial position
        
        %% output state
        dy % position derivative
        dz % scaled velocity derivative
        dx % phase variable derivative

    end
end
