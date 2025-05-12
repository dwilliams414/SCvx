classdef LinearizedModel
    %LinearizedModel: Class for formulating linearized model of nonlinear
    %dynamics under control inputs.  
    
    properties (SetAccess = private)
        f_nl  % Initial dynamics formulation @(t, x, u_k)
        dfdx  % Derivative function of dynamics wrt state @(t, x, u_k)
        dfdu  % Derivative function of dynamics wrt ctrl  @(t, x, u_k)
        n_x   % Number of state variables
        n_u   % Number of control variables
        state_indices
        stm_indices
        B_indices
        c_indices
    end

    methods 
        function obj = LinearizedModel(f_nl, dfdx, dfdu, n_x, n_u)
            %LinearizedModel: Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                f_nl (1, 1) function_handle
                dfdx (1, 1) function_handle
                dfdu (1, 1) function_handle
                n_x  (1, 1) {mustBeInteger}
                n_u  (1, 1) {mustBeInteger}
            end
            obj.f_nl = f_nl;
            obj.dfdx = dfdx;
            obj.dfdu = dfdu;
            obj.n_x = n_x;
            obj.n_u = n_u;

            obj.state_indices   = 1:obj.n_x;
            obj.stm_indices     = obj.n_x+1:obj.n_x+obj.n_x^2;
            obj.B_indices       = obj.stm_indices(end)+1:obj.stm_indices(end)+obj.n_x*obj.n_u;
            obj.c_indices       = obj.B_indices(end)+1:obj.B_indices(end)+obj.n_x;

            % Validate dynamics function
            test_x = randn(n_x, 1);
            test_u = randn(n_u, 1);
            f_test = f_nl(0, test_x, test_u);
            dfdx_test = dfdx(0, test_x, test_u);
            dfdu_test = dfdu(0, test_x, test_u);

            if length(f_test) ~= n_x || ...
                    size(dfdx_test, 1) ~= n_x ||...
                    size(dfdu_test, 1) ~= n_x ||...
                    size(dfdx_test, 2) ~= n_x ||...
                    size(dfdu_test, 2) ~= n_u
                error("Invalid dynamical model!  Dimension mismatch");
            end
        end

        function f = dynamics4discretization(obj)
            arguments (Output)
                f (1, 1) function_handle
            end
            f = @(t, Y, u_k) linearized_with_zoh(t, Y, u_k);
            function dYdt = linearized_with_zoh(t, Y, u_k)
                if (size(u_k) ~= [obj.n_u, 1])
                    error("Invalid control specified!");
                end
                dYdt = zeros(size(Y));

                % Reference State Derivatives
                x = Y(1:obj.n_x);
                dxdt = obj.f_nl(t, x, u_k);
                dYdt(obj.state_indices) = dxdt;

                % STM Derivatives
                stm = reshape(Y(obj.stm_indices), [obj.n_x obj.n_x]);
                A_t = obj.dfdx(t, x, u_k);
                dphidt = A_t*stm;
                dYdt(obj.stm_indices) = dphidt(:);

                % Requirements for B_k
                B_t = obj.dfdu(t, x, u_k);
                B_prod = inv(stm)*B_t;
                dYdt(obj.B_indices) = B_prod(:);

                % Requirements for c_k
                c_t = obj.f_nl(t, x, u_k)-A_t*x-B_t*u_k;
                dYdt(obj.c_indices) = inv(stm)*c_t;
            end
        end
    end
end