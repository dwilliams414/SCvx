classdef SCvxStarFast
    %SCVXSTAR2 Second implementation of the SCvx* framework, to hopefully
    %diagnose where some issues may be occuring in a negative deltaL.
    
    properties

        % ConvexSubproblem to solve at each iteration
        cvxSubproblem (1, 1) CVXFastSubproblem

        % Maximimum trust region radius (rMax)
        rMax (1, 1) double = 0.03;

        % Minimum trust region radius (rMin)
        rMin (1, 1) double = 1e-7;

        % Initial Trust region radius
        initTRRadius (1, 1) double = 0.01;

        % Optimality reduction Factor (gamma)
        gamma (1, 1) double

        % Trust region growth factor (alpha_2)
        alpha2 (1, 1) double {mustBeGreaterThan(alpha2, 1)} = 1.2;

        % Trust region reduction factor (alpha_1)
        alpha1 (1, 1) double {mustBeGreaterThan(alpha1, 1)} = 1.5;

        % Use Fast Solver (TODO)
        useFast (1, 1) {mustBeInteger} = 0;

        % Feasibility Tolerance
        tolFeas (1, 1) double = 5e-9;

        % Optimality Tolerance
        tolOpt (1, 1) double = 1e-6;

        % Max Penalty weight, w
        wMax (1, 1) double

        % Initial Penalty weight
        wInit (1, 1) double

        % Maximum Iterations
        maxIter (1, 1) {mustBeInteger} = 100;

        % Weight increase factor (beta)
        beta (1, 1) double;

        % Minimum Nonlinear Improvement Ratio (deltaJ/deltaL) to Accept step (rho0).  
        rho0 (1, 1) double = 0;

        % Maximum Improvement Ratio to NOT reduce TR Radius
        rho1 (1, 1) double = 0.05;

        % Minimum Improvement Ratio to Increase TR Radius
        rho2 (1, 1) double = 0.4;
    end
    
    methods
        function obj = SCvxStarFast(problem, opts)
            %SCVXSTAR2 Construct an instance of this class
            %   Inputs:
            %       problem: ConvexSubproblem object
            %   Options: TODO
            %   Outputs:
            %       SCvxStar2 object.
            arguments
                problem (1, 1) ConvexSubproblem
                opts.wMax (1, 1) double = 1e6;
                opts.wInit (1, 1) double = 1000;
                opts.weightIncreaseFactor = 1.25
                opts.optimalityReductionFactor (1, 1) double = 0.8;
            end
            obj.cvxSubproblem = problem;
            obj.wMax = opts.wMax;
            obj.wInit = opts.wInit;
            obj.beta = opts.weightIncreaseFactor;
            obj.gamma = opts.optimalityReductionFactor;
        end
        
        function zOpt = solve(obj, zInit)
            %solveStandard: Solve the optimization problem using the SCvx*
            %framework (slow building, no optimizer object in YALMIP)
            %   Inputs:
            %       zInit:      Initial guess for optimization variables
            %   Outputs:
            %       zOpt:       Optimized solution
            arguments
                obj (1, 1) SCvxStarFast
                zInit (:, 1) double
            end
            % Set Current Values for Quantities - Iteration Control:
            zRef = zInit;                   % Current reference opt vars
            currentW = obj.wInit;           % Current penalty weight
            currentTR = obj.initTRRadius;   % Current TR Radius

            optMetric = inf;                % Current Optimality
            feasMetric = inf;               % Current Feasibility
            delta = inf;                    % Current opt improve value

            numIter = 0;                    % Number of iterations

            pRef = obj.cvxSubproblem.build_params(zRef); % Ref param values

            % Get initial g_lin, h_lin to size lagrange multipliers
            glinInit = obj.cvxSubproblem.g_linearized(zRef, pRef);
            hlinInit = obj.cvxSubproblem.h_linear(zRef, pRef);
            currentLambda = zeros(size(glinInit));
            currentMu = zeros(size(hlinInit));

            % Declare SDP Variables for optimization
            zSDP = obj.cvxSubproblem.build_sdpvar(zRef); % Opt vars
            pSDP = sdpvar(length(pRef), 1); % Parameter vector
            zRefSDP = sdpvar(length(zRef), 1); % Reference z values (for TR)
            wSDP = sdpvar(1); % Penalty weight as SDPvar
            lambdaSDP = sdpvar(1); % lambda as sdpvar
            muSDP = sdpvar(length(currentMu, 1)); % mu as sdpvar
            trSDP = sdpvar(1); % trust region as sdpvar
            xi = sdpvar(length(currentLambda), 1); % slack sdpvar
            zeta = sdpvar(length(currentMu), 1); % slack sdpvar


            % Build the Constraints
            gTilde = obj.cvxSubproblem.g_linearized(zSDP, pSDP);
            hTilde = obj.cvxSubproblem.h_linear(zSDP, pSDP);
            gAffine = obj.cvxSubproblem.g_affine(zSDP);
            hConvex = obj.cvxSubproblem.h_cvx(zSDP);

            constraints = [gAffine == 0; hConvex <= 0];

            if ~isempty(gTilde)
                xi = sdpvar(length(gTilde), 1);
                constraints = [constraints; gTilde - xi == 0];
            else
                xi = [];
            end

            if ~isempty(hTilde)
                zeta = sdpvar(length(hTilde), 1);
                constraints = [constraints; hTilde - zeta <= 0];
                constraints = [constraints; zeta >= 0];
            else
                zeta = [];
            end
            
            % Append Trust Region Constraints
            constraints = [constraints; norm(zSDP-zRefSDP, inf) <= trSDP];

            % Define the Objective
            L = obj.augmentedLagrangianSDP(zSDP, xi, zeta, wSDP, ...
                lambdaSDP, muSDP);

            % Solver Settings
            opts = sdpsettings('solver', 'mosek');
            opts.mosek.MSK_DPAR_BASIS_TOL_X = 1e-9;
            opts.mosek.MSK_DPAR_BASIS_TOL_S = 1e-9;
            opts.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS = 1e-11;
            opts.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS = 1e-11;
            opts.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-11;

            % Build Optimizer
            cvxSubproblemOpt = optimizer(constraints, L, opts, ...
                [pSDP; zRefSDP; wSDP; lambdaSDP; muSDP; trSDP], ...
                [zSDP; xi; zeta]);

            sol = cvxSubproblemOpt([pRef; zRef; currentW; currentLambda; currentMu; currentTR]);

        %     while (optMetric > obj.tolOpt || feasMetric > obj.tolFeas) && ...
        %             numIter < obj.maxIter
        % 
        %         numIter = numIter + 1;
        % 
        % 
        % 
        % 
        % 
        % 
        %         % Solve the CVX Subproblem
        %         L = obj.augmentedLagrangianSDP(zSDP, xi, zeta, currentW, ...
        %             currentLambda, currentMu);
        % 
        % 
        % 
        %         sol = optimize(constraints, L, opts);
        % 
        %         % Evaluate Nonlinear Constraints - Current Iteration
        %         gCurrent = obj.cvxSubproblem.g_nonlinear(value(zSDP));
        %         hCurrent = obj.cvxSubproblem.h_nonlinear(value(zSDP));
        %         Jcurrent = obj.augmentedLagrangian(value(zSDP), ...
        %             gCurrent, hCurrent, currentW, currentLambda, ...
        %             currentMu);
        % 
        %         % Evaluate Nonlinear Constraints - Previous Iteration
        %         gPrev = obj.cvxSubproblem.g_nonlinear(zRef);
        %         hPrev = obj.cvxSubproblem.h_nonlinear(zRef);
        %         Jprev = obj.augmentedLagrangian(zRef, gPrev, hPrev, ...
        %             currentW, currentLambda, currentMu);
        % 
        %         % Cost Function Deltas
        %         deltaJCurrent = Jprev-Jcurrent;
        %         deltaLCurrent = Jprev-value(L);
        %         optMetric = abs(deltaJCurrent); % Optimality metric
        % 
        %         if (deltaLCurrent < 0)
        %             warning("This shouldn't be negative...attempt to " + ...
        %                 "diagnose!");
        %         end
        % 
        %         % Evaluate Feasibility
        %         feasVector = [gCurrent; max(0, hCurrent)];
        %         feasMetric = norm(feasVector, inf);
        % 
        %         % Evaluate step improvement
        %         if deltaLCurrent == 0
        %             currentRho = 1;
        %         else
        %             currentRho = deltaJCurrent/deltaLCurrent;
        %         end
        % 
        %         fprintf("Step Results:\n");
        %         fprintf("Iteration Optimality (deltaJ): %.4g\n", deltaJCurrent);
        %         fprintf("Iteration Feasibility: %.4g\n", feasMetric);
        %         fprintf("Iteration deltaL: %.4g\n", deltaLCurrent);
        %         fprintf("Current Trust Region: %.4g\n", currentTR);
        %         fprintf("Current rho: %.4g\n", currentRho);
        %         fprintf("Current w: %.4g\n", currentW);
        %         fprintf("Current Optimality Criteria to Update (delta): %.4g\n", delta);
        %         fprintf("YALMIP Time: %.4g\n", sol.yalmiptime)
        % 
        %         % Assess Step Improvement - keep results
        %         if currentRho >= obj.rho0
        %             zRef = value(zSDP); % Update the reference...good step
        % 
        %             % Assess if we should update weights - we have improved
        %             % optimality enough
        %             if abs(deltaJCurrent) < delta
        %                 fprintf("Updated Multipliers!\n");
        %                 % Update Multipliers
        %                 currentLambda = currentLambda + currentW*gCurrent;
        %                 currentMu = max(0, currentMu + currentW*hCurrent);
        %                 currentW = min(obj.beta*currentW, obj.wMax);
        % 
        %                 % Update Stationary Tolerance
        %                 if isinf(delta)
        %                     delta = abs(deltaJCurrent);
        %                 else
        %                     delta = obj.gamma * delta;
        %                 end
        %             end
        %         end
        % 
        %         % Update Trust Region
        %         if currentRho < obj.rho1
        %             fprintf("Reduced Trust Region\n");
        %             currentTR = max(currentTR/obj.alpha1, obj.rMin);
        %         elseif currentRho > obj.rho2
        %             fprintf("Expanded Trust Region\n");
        %             currentTR = min(obj.alpha2*currentTR, obj.rMax);
        %         end
        %     end
        %     zOpt = value(zSDP);
        % end

        zOpt = solveFast(obj, zInit);

        function L = augmentedLagrangianSDP(obj, zSDP, xi, zeta, w, lambda, mu)
            % AUGMENTEDLAGRANGIANSDP: Augmented lagrangian that for which
            % approximiate solution is obtained at each iteration to result
            % in marching to true optimal of desired problem.  Defined by
            % Equations (4) and (12) in Oguri
            %   Inputs:
            %       obj:    Instance
            %       zSDP:   SDP Variables (YALMIP) for optimization
            %       xi:     SDP Variables (YALMIP) for relaxation of
            %               linearized equality constraints
            %       zeta:   SDP Variables (YALMIP) for relaxation of
            %               convexified inequality constraints
            %       w:      Current penalty weight
            %       lambda: Current lagrange multiplier guess for
            %               linearized equality constraint violation
            %       mu:     Current lagrange multiplier guess for
            %               inequality constraint violation.
            %   Outputs:    Augmented lagrangian function
            f0 = obj.cvxSubproblem.cost_fcn(zSDP);

            L = f0 + dot(lambda, xi) + w/2 * dot(xi, xi) + dot(mu, zeta) + ...
                w/2 * dot(zeta, zeta); % dot(zeta, zeta) because we constraint zeta >= 0!
        end

        function J = augmentedLagrangian(obj, z, gNL, hNL, w, lambda, mu)
            % AUGMENTEDLAGRANGIAN: Evaluate the augmented lagrangian
            % leveraging the nonlinear evaluations of the constraint
            % functions.
            %   Inputs:
            %       obj:        Instance
            %       z:          Value of optimization variables
            %       gNL:        Nonlinear equality constraint eval
            %       hNL:        Non-convex inequality constraint eval
            %       w:          Current penalty weight
            %       lambda:     Current lagrange multiplier estimate
            %       mu:         Current lagrange multiplier estimate (ineq)
            f0 = obj.cvxSubproblem.cost_fcn(z);

            J = f0 + dot(lambda, gNL) + w/2*dot(gNL, gNL) + ...
                dot(mu, hNL) + w/2 * dot(max(0, hNL), max(0, hNL));
        end
    end
end

