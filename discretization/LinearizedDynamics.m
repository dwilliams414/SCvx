classdef LinearizedDynamics
    %LINEARIZEDDYNAMICS: Provide utilities for integrating linearized ODEs
    %as is necessary for discretization with ZOH.  Based upon the method of
    %discretization covered in AAE 590 (ACA) by Prof. Oguri and the
    %MathSpec I wrote up in Overleaf on 05-12-2025.  
    
    properties (SetAccess = protected)
        f_nl
        dfdx
        dfdu
        dfdp
        nx
        nu
        np
        state_indices
        stm_indices
        B_indices
        C_indices
        d_indices
    end
    
    methods
        function obj = LinearizedDynamics(f_nl, dfdx, dfdu, dfdp, ...
                nx, nu, np)
            %LINEARIZEDDYNAMICS Constructor.
            %   Inputs:
            %       f_nl:   Nonlinear dynamics s.t. xdot = f(x, u, p)
            %       dfdx:   Derivative of dynamics wrt state dfdx(x, u, p)
            %       dfdu:   Derivative of dynamics wrt ctrl  dfdu(x, u, p)
            %       dfdp:   Derivative of dynamics wrt param dfdp(x, u, p)
            %       nx:     State dimension 
            %       nu:     Control dimension
            %       np:     Parameter dimension
            %   Outputs:
            %       obj:    Instance of LinearizedDynamics
            arguments
                f_nl (1, 1) function_handle
                dfdx (1, 1) function_handle
                dfdu (1, 1) function_handle
                dfdp (1, 1) function_handle
                nx (1, 1) {mustBeInteger}
                nu (1, 1) {mustBeInteger}
                np (1, 1) {mustBeInteger}
            end

            % Validate Dimensions
            sample_state = zeros(nx, 1);
            sample_ctrl = zeros(nu, 1);
            sample_param = zeros(np, 1);

            A_test = dfdx(sample_state, sample_ctrl, sample_param);
            B_test = dfdu(sample_state, sample_ctrl, sample_param);
            C_test = dfdp(sample_state, sample_ctrl, sample_param);

            if ~all(size(A_test) == [nx nx])
                error("Invalid specification for dfdx!");
            end

            if ~all(size(B_test) == [nx nu])
                error("Invalid specification for dfdu!");
            end

            if ~all(size(C_test) == [nx np])
                error("Invalid specification for dfdp!");
            end

            % If all checks passed, assign values
            obj.f_nl = f_nl;
            obj.dfdx = dfdx;
            obj.dfdu = dfdu;
            obj.dfdp = dfdp;
            obj.nx = nx;
            obj.nu = nu;
            obj.np = np;

            % Get Indices for indexing into augmented state vector
            obj.state_indices = 1:nx;
            obj.stm_indices = obj.state_indices(end)+1:...
                obj.state_indices(end)+obj.nx^2;
            obj.B_indices = obj.stm_indices(end)+1:...
                obj.stm_indices(end)+obj.nx*obj.nu;
            obj.C_indices = obj.B_indices(end)+1:...
                obj.B_indices(end)+obj.nx*obj.np;
            obj.d_indices = obj.C_indices(end)+1:obj.C_indices(end)+obj.nx;
        end

        function f_discretize = dynamics_for_integration(obj, u_ref, p_ref)
            %dynamics_for_integration: Returns the ODE that will allow for
            %propagation to generate the discrete time matrices Ak, Bk, Ck,
            %and dk
            %   Inputs:
            %       obj:    Instance of LinearizedDynamics
            %       uref:   uref(x, p)->control
            %       pref:   Parameter vector np x 1
            %   Ouptus:
            %       f_discretize:   Dynamics function for discretization
            %                       with input signature f_disc(y)
            arguments
                obj (1, 1) LinearizedDynamics
                u_ref (1, 1) function_handle
                p_ref (:, 1) double
            end
            
            % Validate uref function
            sample_state = randn(obj.nx, 1);
            sample_ctrl = randn(obj.nu, 1);

            u_chk = u_ref(sample_state, p_ref);
            if ~all(size(u_chk) == [obj.nu, 1])
                error("Invalid uref function specified!");
            end

            if ~all(size(p_ref) == [obj.np, 1])
                error("Invalid parameter vector specified!");
            end

            % Now build dynamics fcn to return
            function derivative = dynamics(y)
                derivative = zeros(obj.nx+obj.nx^2+obj.nx*obj.nu+...
                    obj.nx*obj.np+obj.nx, 1);

                % Calculate the present control
                x = y(1:obj.nx);
                u = u_ref(x, p_ref);

                % Unpack other variables
                phi = reshape(y(obj.stm_indices), obj.nx, obj.nx);

                % State Derivative
                derivative(obj.state_indices) = obj.f_nl(x, u, p_ref);

                % STM Derivative
                phidot = obj.dfdx(x, u, p_ref)*phi;
                derivative(obj.stm_indices) = phidot(:);

                % B Derivative
                B_deriv = inv(phi)*obj.dfdu(x, u, p_ref);
                derivative(obj.B_indices) = B_deriv(:);

                % C Derivative
                C_deriv = inv(phi)*obj.dfdp(x, u, p_ref);
                derivative(obj.C_indices) = C_deriv(:);

                % d Derivative
                d = obj.f_nl(x, u, p_ref)-obj.dfdx(x, u, p_ref)*x-...
                    obj.dfdu(x, u, p_ref)*u-obj.dfdp(x, u, p_ref)*p_ref;
                d_deriv = inv(phi)*d;

                derivative(obj.d_indices) = d_deriv;
            end

            % Dynamics fcn
            f_discretize = @(t, y) dynamics(y);
        end

        function [Ak, Bk, Ck, dk, xkp1] = integrate_discretized(obj, xk, ...
                u_ref, p_ref, tspan, opts)
            %integrate_discretized: Integrate initial state, xk, for t in
            %tk to tk+1 to generate the associated discretization matrices
            %and final state
            %   Inputs:
            %       xk:     nx X 1 initial state
            %       u_ref:  u_ref(x, p_ref)-> u (control)
            %       p_ref:  Reference parameter vector
            %       tspan:  Timespan of integration, tk->tk+1
            %  Options:
            %       integrator:     Choice of integrator to use
            %       abs_tol:        Absolute tolerance for integration
            %       rel_tol:        Relative tolerance for integration
        arguments
            obj (1, 1) LinearizedDynamics
            xk (:, 1) double
            u_ref (1, 1) function_handle
            p_ref (:, 1) double
            tspan (1, :) double
            opts.integrator (1, 1) function_handle = @ode89;
            opts.abs_tol (1, 1) double = 1e-12;
            opts.rel_tol (1, 1) double = 1e-12;
        end
            f_dynamics = obj.dynamics_for_integration(u_ref, p_ref);
            odeopts = odeset('RelTol', opts.rel_tol, 'AbsTol', ...
                opts.abs_tol);


            phi0 = eye(obj.nx);
            B_init = zeros(obj.nx, obj.nu);
            C_init = zeros(obj.nx, obj.np);
            d_init = zeros(obj.nx, 1);
            y0 = [xk; phi0(:); B_init(:); C_init(:); d_init];

            [t, x] = opts.integrator(f_dynamics, tspan, y0, odeopts);

            Ak = reshape(x(end, obj.stm_indices), obj.nx, obj.nx);
            Bk = Ak * reshape(x(end, obj.B_indices), obj.nx, obj.nu);
            Ck = Ak * reshape(x(end, obj.C_indices), obj.nx, obj.np);
            dk = Ak * reshape(x(end, obj.d_indices), obj.nx, 1);
            xkp1 = x(end, 1:obj.nx)';
        end
    end
end

