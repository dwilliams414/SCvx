classdef (Abstract) CVXFastSubproblem
    %CVXFastSubproblem: Specifies necessary components for constructing the
    %convex subproblem to solve at each iteration of the SCvx(*) algorithm
    %using an optimizer object in YALMIP.  In this case, a parameter
    %generation function is required that is a function of z_ref, but all
    %evaluations of g_lin, etc., must be a function of p, as opposed to
    %z_ref explicitly.

    properties (SetAccess = public)
    end

    methods (Abstract)
        % g_linearized: Evaluate linearized equality constraints for a
        % given convex subproblem.
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %       z_ref:  Reference optimization variable vector for
        %               linearization
        %   Outputs:
        %       Evaluate g_tilde(z) for SCvx* framework for the
        %       linearization specified by z_ref
        g_lin = g_linearized(obj, z, pVec);

        % g_nonlinear:  Evaluate the nonlinear equality constraints for a
        % specified optimization variable vector, z.  Used for evaluating,
        % e.g., nonlinear improvement/feasibility in SCvx* framework.
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %   Outputs:
        %       g_nl(z)
        g_nl = g_nonlinear(obj, z);

        % g_affine:  Evaluate the affine equality constraints for a
        % specified optimization variable vector, z.  Used for evaluating,
        % e.g., nonlinear improvement/feasibility in SCvx* framework.
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %   Outputs:
        %       g_nl(z)
        g_aff = g_affine(obj, z);

        % h_linear: Evaluate the linearized inequality constraints about a
        % reference set of optimization variables, z_ref
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %       z_ref:  Reference optimization variable vector for
        %               linearization
        %   Outputs:
        %       Evaluate h_tilde(z) for SCvx* framework for the
        %       linearization specified by z_ref
        h_lin = h_linear(obj, z, pVec);


        % h_nonlinear:  Evaluate the nonlinear inequality constraints for a
        % specified optimization variable vector, z.  Used for evaluating,
        % e.g., nonlinear improvement/feasibility in SCvx* framework.
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %   Outputs:
        %       h_nonlinear(z)
        h_nl = h_nonlinear(obj, z);

        % h_cvx:  Evaluate the nonlinear inequality constraints for a
        % specified optimization variable vector, z.  Used for evaluating,
        % e.g., nonlinear improvement/feasibility in SCvx* framework.
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %   Outputs:
        %       h_cvx(z)
        h_convex = h_cvx(obj, z);

        % cost_fcn: Evaluate the cost function to minimize, f_0(z) in the
        % SCvx* notation.
        %   Inputs:
        %       z:      Optimization variable vector (n x 1)
        %   Outputs:
        %       f0_z:   Cost function value (1, 1)
        f0 = cost_fcn(obj, z);

        % build_sdpvar: Build the YALMIP SDPVAR that represents the desired
        % subproblem.  We do this here because we want to be able to
        % incorporate matrices and LinearMatrixInequality (LMI) constraints
        %   Inputs:
        %       obj:        Instance of ConvexSubproblem
        %       z_ref:      Reference optimization variable values (if
        %                   needed)
        %   Outputs:
        %       z_sdp:          (Vector) SDPVAR
        z_sdp = build_sdpvar(obj, z_ref);

        p = build_params(obj, z_ref);
    end

    methods
        function g_all_nonlinear = g_all_nl(obj, z)
            % Evaluate all equality constraints (nonlinear)
            g_all_nonlinear = [obj.g_nonlinear(z); obj.g_affine(z)];
        end

        function g_all_linear = g_all_lin(obj, z, z_ref)
            % Evaluate all equality constraints (linear)
            g_all_linear = [obj.g_linearized(z, z_ref); obj.g_affine(z)];
        end

        function h_all_nonlinear = h_all_nl(obj, z)
            % Evaluate all inequality constraints (nonlinear)
            h_all_nonlinear = [obj.h_nonlinear(z); obj.h_cvx(z)];
        end

        function h_all_nonlinear = h_all_lin(obj, z, z_ref)
            % Evaluate all inequality constraints (linear)
            h_all_nonlinear = [obj.h_linear(z, z_ref); obj.h_cvx(z)];
        end
    end
end