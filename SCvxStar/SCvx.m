classdef SCvx
    %SCVX: Implementation of the SCvx(*) framework for successive convex
    %optimization with a feasibility guarantee. For mathematical
    %specification and definition of quantities, refer to Oguri 2023
    %available here: https://arxiv.org/pdf/2304.14564
    
    properties
        Property1
    end
    
    methods
        function obj = SCvx(inputArg1,inputArg2)
            %SCVX Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

