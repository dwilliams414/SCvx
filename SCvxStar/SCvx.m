classdef SCvx
    %SCVX: Implementation of the SCvx(*) framework for successive convex
    %optimization with a feasibility guarantee. For mathematical
    %specification and definition of quantities, refer to Oguri 2023
    %available here: https://arxiv.org/pdf/2304.14564
    
    properties (SetAccess = protected)
        % Convex subproblem to solve at each iteration of algorithm.
        % Instance of the ConvexSubproblem object
        cvx_subproblem

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
            opts.tol_opt = tol_opt;
        end
        
        function outputArg = solve(obj, z_ref)
            %solve: Solve the optimization problem using the SCvx*
            %framework
            %   Inputs:
            %       obj:    Self
            %       z_ref:  Reference optimization variables for first
            %               iteration
            %   Outputs:
            %       (TBD): Some kind of solution structure
            
        end
    end
end

