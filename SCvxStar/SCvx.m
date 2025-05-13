classdef SCvx
    %SCVX: Implementation of the SCvx(*) framework for successive convex
    %optimization with a feasibility guarantee. For mathematical
    %specification and definition of quantities, refer to Oguri 2023
    %available here: https://arxiv.org/pdf/2304.14564
    
    properties (SetAccess = protected)
        % Convex subproblem to solve at each iteration of algorithm.
        % Instance of the ConvexSubproblem object
        cvx_subproblem (1, 1) % ConvexSubproblem % For coding rn, remove

        % Equality constraint Lagrange multipliers for augmented lagrangian
        lambda (:, 1) double

        % Inequality constraint Lagrange multipliers for augmented
        % lagrangian
        lambda_ineq (:, 1) double

        % Current penalty weight, w
        w_current (1, 1) double
  
    end

    properties (SetAccess = public)
        % Feasibility tolerance for acceptable solution (stopping criteria)
        tol_feas (1, 1) double = 1e-9;

        % Optimality tolerance for acceptable solution (stopping criteria)
        tol_opt (1, 1) double = 1e-6;

        % Initial penalty weight (w in Oguri 2023)
        w_init(1, 1) double = 100;

        % Growth factor for penalty weight, w, denoted beta in SCvx*
        beta_factor (1, 1) double = 4;

        % Reduction factor for delta, determining when multipliers are
        % updated.  See Eqn 16 in Oguri.
        gamma = 0.75;

        % Minimum trust region radius
        min_tr_radius = 1e-9;

        % Maximum trust region radius
        max_tr_radius = 0.1;

        % Trust region reduction factor (> 1)
        alpha1 = 1.5;

        % Trust region growth factor (> 1)
        alpha2 = 1.2;

        % Initial trust region radius
        init_tr_radius = 0.01;

        % Maximum number of iterations
        max_iter = 100;
    end
    
    methods
        function obj = SCvx(problem, opts)
            %SCVX Construct an instance of this class
            %   Inputs:
            %       problem:    ConvexSubproblem to solve.  Note that this
            %                   convex subproblem does NOT include the
            %                   trust region constraints/augmented
            %                   Lagarangian that are to be managed by SCvx
            %   Options: (TODO)

            %   Outputs:
            %       obj:        SCvx object
            arguments
                problem (1, 1) ConvexSubproblem
                opts.tol_feas (1, 1) double = 1e-9;
                opts.tol_opt (1, 1) double = 1e-6;
            end

            obj.cvx_subproblem = problem;
            obj.tol_feas = opts.tol_feas;
            opts.tol_opt = opts.tol_opt;
            
        end
        
        function z_opt = solve(obj, z_ref)
            %solve: Solve the optimization problem using the SCvx*
            %framework
            %   Inputs:
            %       obj:    Self
            %       z_ref:  Reference optimization variables for first
            %               iteration
            %   Outputs:
            %       z_opt: Optimized value at final iteration
            fprintf("Time to implement this!!");

            % Initial weight setup
            obj.w_current = obj.w_init;

            % Stopping Conditions
            feasibility_metric = inf;
            deltaJ = inf;
            num_iter = 0;

            while (feasibility_metric > obj.tol_feas...
                    || deltaJ > obj.tol_opt) && num_iter < obj.max_iter

                % Construct SDP variables
                z_sdp = obj.cvx_subproblem.build_sdpvar(z_ref);

                % Build the LHS of the Equality Constraints
                g_lin = obj.cvx_subproblem.g_linearized(z_sdp, z_ref);
                g_affine = obj.cvx_subproblem.g_affine(z_sdp);

                % Build the LHS of the Inequality Constraints
                h_lin = obj.cvx_subproblem.h_linear(z_sdp, z_ref);
                h_cvx = obj.cvx_subproblem.h_cvx(z_sdp);

                % Construct overall constraint vector for convex iteration
                xi = sdpvar(length(g_lin), 1);
                zeta = sdpvar(length(h_lin), 1);

                constraints = [g_lin == xi; ...
                    h_lin <= zeta; zeta >= 0; g_affine == 0; h_cvx <= 0];

                % Build Augmented Lagrangian
                if (num_iter == 0)
                    obj.lambda = zeros(length(g_lin)+length(g_affine), 1);
                    obj.lambda_ineq = zeros(length(h_lin)+length(h_cvx), 1);
                end
                L = obj.build_augmented_lagrangian(z_sdp, ...
                    [g_lin; g_affine], [h_lin; h_cvx]);

                % Solve Iteration
                sol = optimize(constraints, L);

                fprintf("Iteration Done!\n");
                num_iter = num_iter + 1;
            end
        end

        function aug_lagrangian = build_augmented_lagrangian(obj, z,...
                g, h)
            arguments (Input)
                obj (1, 1) SCvx
                z (:, 1)
                g (:, 1)
                h (:, 1)
            end
            arguments (Output)
                aug_lagrangian (1, 1)
            end
            f0 = obj.cvx_subproblem.cost_fcn(z);

            h_plus = max(0, h);
            aug_lagrangian = f0 + obj.lambda'*g + ...
                obj.lambda_ineq'*h + ...
                obj.w_current/2 * (g'*g + h_plus'*h_plus);
        end
    end
end

