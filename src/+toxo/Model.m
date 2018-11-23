classdef Model
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        order
        pt_s
    end
    
    methods (Access = protected, Static = true)
        
        
        function [c, ceq] = hconstraint(pt, h)
            c = [];
            ceq = pt.heritability() - h;
        end
    end
    
    methods
        function obj = Model(path)
            %MODEL Construct an instance of this class
            %   Detailed explanation goes here
            
            fid = fopen(path, 'r');
            content = textscan(fid, "%s %s");
            fclose(fid);
            [~, index] = sort(content{1});
            obj.pt_s = str2sym(content{2}(index));
            obj.order = length(content{1}{1}) / 2;
        end
        
        function p = genotype_probabilities(obj, maf)
            p = prod(toxo.nfold([(1 - maf)^2, 2 * maf * (1 - maf), maf^2], obj.order), 2);
        end
        
        function pt = find_max_prevalence(obj, maf, h)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % Create minimization problem
            
            highest_degree_monomial = obj.pt_s(end);
            bfunction = matlabFunction( ...
                rhs(isolate(highest_degree_monomial == 1, sym('b'))) ...
            );
            lb = 0;
            ub = 1;
            x0 = 0.01;
            nonlcon = @(x) toxo.Model.hconstraint( ...
                toxo.PT(obj, maf, x, bfunction(x)), h ...
            );
            options = optimoptions('fmincon', ...
                'Display', 'none', ...
                'StepTolerance', 1e-120, ...
                'ConstraintTolerance', 1e-3 ...
            );
            [alpha, beta, exitflag, output] = fmincon( ...
                bfunction, x0, [], [], [], [], lb, ub, nonlcon, options ...
            );
            if exitflag < 0
                ME = MException("toxo.Model:no_solution", output.message);
                throw(ME);
            end
            
            % Obtain the resulting heritability and penetrance table
            pt = toxo.PT(obj, maf, alpha, beta);
        end
    end
end

