classdef SCvxStar2
    %SCVXSTAR2 Second implementation of the SCvx* framework, to hopefully
    %diagnose where some issues may be occuring in a negative deltaL.
    
    properties

        % ConvexSubproblem to solve at each iteration
        cvxSubproblem (1, 1)

        % Maximimum trust region radius (rMax)
        rMax (1, 1) double = 0.1;

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
        function obj = SCvxStar2(problem, opts)
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
        
        function zOpt = solveStandard(obj, zInit)
            %solveStandard: Solve the optimization problem using the SCvx*
            %framework (slow building, no optimizer object in YALMIP)
            %   Inputs:
            %       zInit:      Initial guess for optimization variables
            %   Outputs:
            %       zOpt:       Optimized solution
            arguments
                obj (1, 1) SCvxStar2
                zInit (:, 1) double
            end
            % Set Current Values for Quantities:
            zRef = zInit;
            currentW = obj.wInit;
            currentTR = obj.initTRRadius;
            currentLambda = 0;
            currentMu = 0;

            optMetric = inf;
            feasMetric = inf;
            delta = inf; % delta from SCVx*

            numIter = 0;
            while optMetric > obj.tolOpt && feasMetric > obj.tolFeas && ...
                    numIter < obj.maxIter

                numIter = numIter + 1;
                yalmip('clear'); % Clear existing SDP variables

                % Build the optimization variable vector
                zSDP = obj.cvxSubproblem.build_sdpvar(zRef);

                % Generate the linearized constraints and convexified
                % inequality constraints
                gTilde = obj.cvxSubproblem.g_linearized(zSDP, zRef);
                hTilde = obj.cvxSubproblem.h_linear(zSDP, zRef);

                % Generate the Affine and Convex constraints
                gAffine = obj.cvxSubproblem.g_affine(zSDP);
                hConvex = obj.cvxSubproblem.h_cvx(zSDP);

                % Initialize the (Subproblem) constraint vector
                constraints = [gAffine == 0; hConvex <= 0];

                if ~isempty(gTilde)
                    xi = sdpvar(length(gTilde, 1));
                    constraints = [constraints; gTilde - xi == 0];
                else
                    xi = [];
                end

                if ~isempty(hTilde)
                    zeta = sdpvar(length(hTilde), 1);
                    constraints = [constraints; hTilde - zeta <= 0];
                    constraints = [constraings; zeta >= 0];
                else
                    zeta = [];
                end

                % Append Trust Region Constraints
                constraints = [constraints; norm(zSDP-zRef, inf) <= currentTR];

                % Solve the CVX Subproblem
                L = obj.augmentedLagrangianSDP(zSDP, xi, zeta, currentW, ...
                    currentLambda, currentMu);

                opts = sdpsettings('solver', 'mosek');
                sol = optimize(constraints, L, opts);

                % Evaluate Nonlinear Constraints - Current Iteration
                gCurrent = obj.cvxSubproblem.g_nonlinear(value(zSDP));
                hCurrent = obj.cvxSubproblem.h_nonlinear(value(zSDP));
                Jcurrent = obj.augmentedLagrangian(value(zSDP), ...
                    gCurrent, hCurrent, currentW, currentLambda, ...
                    currentMu);

                % Evaluate Nonlinear Constraints - Previous Iteration
                gPrev = obj.cvxSubproblem.g_nonlinear(zRef);
                hPrev = obj.cvxSubproblem.h_nonlinear(zRef);
                Jprev = obj.augmentedLagrangian(zRef, gPrev, hPrev, ...
                    currentW, currentLambda, currentMu);
                
                % Cost Function Deltas
                deltaJCurrent = Jprev-Jcurrent;
                deltaLCurrent = Jprev-value(L);
                optMetric = deltaJCurrent; % Optimality metric

                if (deltaLCurrent < 0)
                    warning("This shouldn't be negative...attempt to " + ...
                        "diagnose!");
                end
                
                % Evaluate Feasibility
                feasVector = [gCurrent; max(0, hCurrent)];
                feasMetric = norm(feasVector);

                % Evaluate step improvement
                if deltaLCurrent == 0
                    currentRho = 1;
                else
                    currentRho = deltaJCurrent/deltaLCurrent;
                end

                fprintf("Step Results:\n");
                fprintf("Iteration Optimality (deltaJ): %.4g\n", deltaJCurrent);
                fprintf("Iteration Feasibility: %.4g", feasMetric);
                fprintf("Iteration deltaL: %.4\n", deltaLCurrent);
                fprintf("Current Trust Region: %.4g", currentTR);
                fprintf("Current rho: %.4g\n", currentRho);
                fprintf("Current w: %.4g\n", currentW);
                fprintf("Current Optimality Criteria to Update (delta): %.4g\n", delta);

                % Assess Step Improvement - keep results
                if currentRho >= obj.rho0
                    zRef = value(zSDP); % Update the reference...good step

                    % Assess if we should update weights - we have improved
                    % optimality enough
                    if abs(deltaJCurrent) < delta
                        fprintf("Updated Multipliers!\n");
                        % Update Multipliers
                        currentLambda = currentLambda + currentW*gCurrent;
                        currentMu = max(0, currentMu + currentW*hCurrent);
                        currentW = min(obj.beta*currentW, obj.wMax);

                        % Update Stationary Tolerance
                        delta = obj.gamma * delta;
                    end
                end

                % Update Trust Region
                if currentRho < obj.rho1
                    fprintf("Reduced Trust Region");
                    currentTR = max(currentTR/obj.alpha1, obj.rMin);
                elseif currentRho > obj.rho2
                    fprintf("Expanded Trust Region");
                    currentTR = min(obj.alpha2*currentTR, obj.rMax);
                end
            end
            zOpt = value(zSDP);
        end

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

            L = f0 + dot(lambda, xi) + w/2 * dot(xi, xi) + mu * zeta + ...
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

