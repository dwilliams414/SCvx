classdef SCvx
    %SCVX: Implementation of the SCvx(*) framework for successive convex
    %optimization with a feasibility guarantee. For mathematical
    %specification and definition of quantities, refer to Oguri 2023
    %available here: https://arxiv.org/pdf/2304.14564
    
    properties (SetAccess = protected)
        % Convex subproblem to solve at each iteration of algorithm.
        % Instance of the ConvexSubproblem object
        cvx_subproblem (1, 1) %  ConvexSubproblem % For coding rn, remove

        % Equality constraint Lagrange multipliers for augmented lagrangian
        lambda (:, 1) double

        % Inequality constraint Lagrange multipliers for augmented
        % lagrangian
        lambda_ineq (:, 1) double

        % Current penalty weight, w
        w_current (1, 1) double

        % Current Trust Region Radius
        current_tr_radius (1, 1) double

        % Current delta (see Oguri for definition.  Determines criteria for
        % updating Lagrange multipliers)
        delta_current (1, 1) double

        % Current feasibility metric
        feasibility_metric (1, 1) double

        % Optimality Metric
        optimality_metric (1, 1) double

        % Current iteration number
        num_iter (1, 1) double;
  
    end

    properties (SetAccess = public)
        % Feasibility tolerance for acceptable solution (stopping criteria)
        tol_feas (1, 1) double = 1e-9;

        % Optimality tolerance for acceptable solution (stopping criteria)
        tol_opt (1, 1) double = 1e-6;

        % Initial penalty weight (w in Oguri 2023)
        w_init(1, 1) double = 10;

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
        max_iter = 200;

        % Minimum value for deltaJ/deltaL such that step is accepted
        min_rho_to_accept = 0;

        % If abs(deltaJ/deltaL) > shrink_tr_criteria -> shrink trust region
        shrink_tr_criteria = 0.2;

        % If abs(deltaJ/deltaL) < grow_tr_criteria -> grow trust region
        grow_tr_criteria = 0.1;

        % Maximum allowable penalty weight, w
        max_w = inf;
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
            % 
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

            % Initialize SCvx* Parameters
            obj.w_current = obj.w_init;
            obj.current_tr_radius = obj.init_tr_radius;
            obj.delta_current = inf;

            % Stopping Conditions
            obj.feasibility_metric = inf;
            obj.optimality_metric = inf;
            obj.num_iter = 0;

            while (obj.feasibility_metric > obj.tol_feas...
                    || obj.optimality_metric > obj.tol_opt) && ...
                    obj.num_iter < obj.max_iter

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

                % Append Trust-Region Constraints
                constraints = [constraints; norm(z_sdp-z_ref, 2) ...
                    <= obj.current_tr_radius];

                % Build Augmented Lagrangian
                if (obj.num_iter == 0)
                    obj.lambda = zeros(length(g_lin), 1);
                    obj.lambda_ineq = zeros(length(h_lin), 1);
                end
                
                L = obj.augmented_lagrangian_sdp(z_sdp, xi, zeta);

                opts = sdpsettings('solver', 'mosek');
                % Solve Iteration Problem
                sol = optimize(constraints, L, opts);

                % Compute Cost for Previous Iteration
                g_prev = obj.cvx_subproblem.g_nonlinear(z_ref);
                h_prev = obj.cvx_subproblem.h_nonlinear(z_ref);

                J_prev = obj.augmented_lagrangian(z_ref, g_prev, h_prev);

                % Compute Cost for Current Iteration
                g_current = obj.cvx_subproblem.g_nonlinear(value(z_sdp));
                h_current = obj.cvx_subproblem.h_nonlinear(value(z_sdp));

                J_current = obj.augmented_lagrangian(value(z_sdp), ...
                    g_current, h_current);

                % Compute Actual ans Expected Cost Change
                deltaJ = J_prev-J_current;
                deltaL = J_prev-value(L);

                if deltaL < 0
                    warning("This is odd...");
                end

                % Compute Feasibility at Current Iteration
                current_feasibility = norm([g_current; max(0, h_current)], 2);
                norm(value(xi), 2)

                % Calculate rho: observed change/expected
                if deltaL ~= 0
                    rho = deltaJ/deltaL;
                else
                    rho = 1;
                end

                % If we saw improvement in the nonlinear cost, 
                % accept the step
                if rho > obj.min_rho_to_accept
                    % Update Reference
                    z_ref = value(z_sdp);
                    obj.feasibility_metric = current_feasibility;
                    obj.optimality_metric = abs(deltaJ);

                    % Initialize delta_current, if needed
                    if obj.num_iter == 0
                        obj.delta_current = abs(deltaJ);
                    end

                    % Check if we update the Lagrange Multipliers and
                    % penalty weight. Do so if required
                    if abs(deltaJ) <= obj.delta_current
                        fprintf("Updating Multipliers...\n");
                        obj.lambda = obj.lambda + obj.w_current * g_current;
                        obj.lambda_ineq = max(0, obj.lambda_ineq + ...
                            obj.w_current * h_current);
                        obj.w_current = min(obj.max_w,...
                            obj.beta_factor*obj.w_current);

                        obj.delta_current = obj.gamma * obj.delta_current;
                    end
                end

                % Resize Trust Region
                if abs(rho-1) < obj.grow_tr_criteria
                    obj.current_tr_radius = min(obj.max_tr_radius, ...
                        obj.alpha2*obj.current_tr_radius);
                end
                if abs(rho-1) > obj.shrink_tr_criteria
                    obj.current_tr_radius = max(obj.min_tr_radius, ...
                        obj.current_tr_radius/obj.alpha1);
                end
                fprintf("Current Feasibility: %.4g\n", obj.feasibility_metric);
                fprintf("Current Optimality: %.4g\n", obj.optimality_metric);
                fprintf("w: %.4g", obj.w_current);
   
                % Increment Counter
                obj.num_iter = obj.num_iter + 1;
            end

            if obj.num_iter < obj.max_iter
                z_opt = value(z_ref);
            else
                error("Failed to converge!");
            end
        end

        function aug_lagrangian = augmented_lagrangian_sdp(obj, z,...
                xi, zeta)
            % build_augmented_lagrangian: Build the augmented Lagrangian
            % required for SCvx*, following Oguri (2023).  Note that g and
            % h are the FULL equality constraint and inequality constraint
            % violations
            %   Inputs:
            %       obj:    Class instance
            %       z:      Optimization variable vector
            %       g:      Equality constraints (all)
            %       h:      Inequality constraints (all)
            arguments (Input)
                obj (1, 1) SCvx
                z (:, 1)
                xi (:, 1)
                zeta (:, 1)
            end
            arguments (Output)
                aug_lagrangian (1, 1)
            end
            f0 = obj.cvx_subproblem.cost_fcn(z);

            aug_lagrangian = f0 + obj.lambda'*xi + ...
                obj.lambda_ineq'*zeta + ...
                obj.w_current/2 * dot(xi, xi) + obj.w_current/2 * ...
                dot(zeta, zeta);
        end

        function aug_lagrangian = augmented_lagrangian(obj, z,...
                g_nl, h_nl)
            % build_augmented_lagrangian: Build the augmented Lagrangian
            % required for SCvx*, following Oguri (2023).  Note that g and
            % h are the FULL equality constraint and inequality constraint
            % violations
            %   Inputs:
            %       obj:    Class instance
            %       z:      Optimization variable vector
            %       g:      Equality constraints (all)
            %       h:      Inequality constraints (all)
            arguments (Input)
                obj (1, 1) SCvx
                z (:, 1)
                g_nl (:, 1)
                h_nl (:, 1)
            end
            arguments (Output)
                aug_lagrangian (1, 1)
            end
            f0 = obj.cvx_subproblem.cost_fcn(z);

            h_nl_plus = max(0, h_nl);
            aug_lagrangian = f0 + obj.lambda'*g_nl + ...
                obj.lambda_ineq'*h_nl + ...
                obj.w_current/2 * dot(g_nl, g_nl) + obj.w_current/2 * ...
                dot(h_nl_plus, h_nl_plus);
        end
    end
end

