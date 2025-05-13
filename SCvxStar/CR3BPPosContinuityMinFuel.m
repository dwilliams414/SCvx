classdef CR3BPPosContinuityMinFuel < ConvexSubproblem
    %CR3BPPosContinuityMinFuel: Example ConvexSubproblem implementation,
    %solving for a position continuous, minimum fuel trajectory with path
    %constraints.

    properties (SetAccess = private)
        x_ref_init
        nd_times
        stm_propagator % stm_propagator(x0)->[xf, stm_f]
        
        num_nodes; % Number of reference states
        n_x = 6; % State dimension
        n_u = 3; % Control dimension

        state_indices;      % Indices in z (opt variables) corr. to states
        ctrl_indices;    % Indices in z corr. to control inputs
    end

    methods
        function obj = CR3BPPosContinuityMinFuel(x_ref_init, nd_times, ...
                mass_ratio)
            %EXAMPLESUBPROBLEM: Define an optimization problem for position
            %continuity while minimizing total impulsive maneuver delta-V.
            %   Inputs:
            %       x_ref_init: (6) X (num_states) array of
            %                   reference states.
            %       nd_times:   Epochs correspondinging to x_ref_init.
            %                   Note that these are not  
            %       mass_ratio: CR3BP mass ratio
            %   Outputs:
            %       obj:    Instance of class
            arguments
                x_ref_init (6, :) double
                nd_times (1, :) double
                mass_ratio (1, 1) double
            end

            if length(nd_times) ~= size(x_ref_init, 2)
                error("Invalid time/state pairing!");
            end

            obj.x_ref_init = x_ref_init;
            obj.nd_times = nd_times;
            
            obj.stm_propagator = CR3BPPosContinuityMinFuel....
                .build_prop_fcn(mass_ratio);

            obj.num_nodes = size(x_ref_init, 2);
            obj.state_indices = 1:obj.num_nodes*obj.n_x;
            obj.ctrl_indices = obj.state_indices(end)+1:...
                obj.state_indices(end)+obj.n_u*(obj.num_nodes-1);
        end

        function g_lin = g_linearized(obj, z, z_ref)
            arguments
                obj (1, 1) CR3BPPosContinuityMinFuel
                z (:, 1)
                z_ref (:, 1)
            end
            A = zeros(obj.n_x*(obj.num_nodes-1));
            B = zeros(size(A, 1), obj.n_u*(obj.num_nodes-1));
            d_vec = zeros(obj.n_x*(obj.num_nodes-1), 1);
            Xkp1 = z(obj.state_indices(7:end));
            Xk = z(obj.state_indices(1:end-6));
            Uk = z(obj.ctrl_indices);

            for k = 1:obj.num_nodes-1
                x0k = z_ref(obj.state_indices(obj.n_x*(k-1)+1:obj.n_x*k));
                uk = z_ref(obj.ctrl_indices(obj.n_u*(k-1)+1:obj.n_u*k));
                tspank = obj.nd_times(k:k+1);

                [xkp1, stmf] = obj.stm_propagator(...
                    [x0k+[zeros(3);eye(3)]*uk], tspank);

                Ak = stmf;
                Bk = Ak*[zeros(3); eye(3)];
                
                A(obj.n_x*(k-1)+1:obj.n_x*k, obj.n_x*(k-1)+1:obj.n_x*k)=Ak;
                B(obj.n_x*(k-1)+1:obj.n_x*k, obj.n_u*(k-1)+1:obj.n_u*k)=Bk;

                d_vec(obj.n_x*(k-1)+1:obj.n_x*k) = xkp1-Ak*x0k-Bk*uk;
    
            end
            g_lin = Xkp1-(A*Xk)-B*Uk-d_vec;
        end

        function g_nl = g_nonlinear(obj, z)
            g_nl = zeros(obj.state_indices(end)-6, 1);
            Xkp1 = z(obj.state_indices(7:end));
            for k = 1:obj.num_nodes-1
                current_state_indices = obj.n_x*(k-1)+1:obj.n_x*k;
                current_ctrl_indices = obj.n_u*(k-1)+1:obj.n_u*k;

                xk = z(obj.state_indices(current_state_indices));
                uk = z(obj.ctrl_indices(current_ctrl_indices));
                tspank = obj.nd_times(k:k+1);

                [xf, ~] = obj.stm_propagator([xk+[zeros(3);eye(3)]*uk], ...
                    tspank);

                g_nl(current_state_indices) = Xkp1(current_state_indices)-...
                    xf;
            end
        end

        function g_aff = g_affine(obj, z)
            g_aff = [];
        end

        function h_convex = h_cvx(obj, z)
            h_convex = [];
            for k = 1:obj.num_nodes
                state_current = ...
                    z(obj.state_indices(obj.n_x*(k-1)+1:obj.n_x*k));

                ref_state = obj.x_ref_init(:, k);
                h_convex = [h_convex; norm(state_current(1:3)-ref_state(1:3), 2)-0.05];
            end
        end

        function h_nl = h_nonlinear(obj, z)
            h_nl = [];
        end

        function h_lin = h_linear(obj, z, zref)
            h_lin = [];
        end

        function f0 = cost_fcn(obj, z)
            u_all = reshape(z(obj.ctrl_indices), 3, []);
            f0 = norm(u_all(:, 1));
            for k = 2:obj.num_nodes-1
                f0 = f0 + norm(u_all(:, k));
            end
        end

        function z_sdp = build_sdpvar(obj, z_ref)
            arguments
                obj (1, 1) CR3BPPosContinuityMinFuel
                z_ref (:, 1) double
            end
            len_z = obj.n_x*obj.num_nodes+obj.n_u*(obj.num_nodes-1);
            z_sdp = sdpvar(len_z, 1);            
        end

    end

    methods (Static)
        function prop_fcn = build_prop_fcn(mass_parameter)
            function [xf, stmf] = propagator(x0, tspan)
                phi0 = eye(6);
                [~, x] = IntegrateCR3BP([x0; phi0(:)], tspan, mass_parameter);
                xf = x(end, 1:6)';
                stmf = state2phimats(x(end, :)');
            end
            prop_fcn = @(x0, tspan) propagator(x0, tspan);
        end
    end
end

