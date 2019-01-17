classdef Model
    %MODEL Abstract representation of a penetrance model as a function
    % of two variables. This class provides methods to obtain numeric
    % penetrance tables derived from the parametric representation given as
    % input.
    
    properties
        name          % Name of the model.
        order         % Number of loci involved in the epistatic model.
        penetrances   % Array of symbolic expressions, representing the epistatic model.
        variables     % List of all variables contained in all symbolic expressions
    end
    
    methods (Access = private, Static = true)
        function [c, ceq] = hconstraint(pt, h)
            c = [];
            ceq = pt.heritability() - h;
        end
        
        function [c,ceq] = pconstraint(pt, p)
            c = [];
            ceq = pt.prevalence() - p;
        end
    end
    
    methods
        function obj = Model(path)
            %MODEL Construct an instance of this class from the given model.
            % M = MODEL(P) reads the model from its text representation in
            % file path P.
            %
            % The input model must be formatted as a plain CSV, with each
            % line of the file P corresponding to a row of the model. The 
            % rows are made of the genotype definition and the probability
            % associated with the given genotype, separated by a coma.
            % Probability is expressed as a function of two variables.
            % Empty lines, as well as lines starting with '#' will be
            % ignored.
            
            [~, obj.name, ~] = fileparts(path);
            fid = fopen(path, 'r');
            content = textscan(fid, "%[^,], %[^,\n]", "CommentStyle", "#", "TextType", "string", "CollectOutput", true, "ReturnOnError", false);
            fclose(fid);
            content = content{1};
            [~, index] = sort(content(:,1));
            obj.penetrances = str2sym(content(index,2));
            obj.variables = symvar(obj.penetrances);
            obj.order = strlength(content(1,1)) / 2;
        end
        
        function pt = find_max_prevalence(obj, maf, h)
            %FIND_MAX_PREVALENCE Compute the penetrance table of the model with the maximum admissible prevalence given its MAF and heritability.
            % PT = FIND_MAX_PREVALENCE(MAF, H) returns a PT instance with the
            % maximum admissible prevalence given the MAF and heritability
            % restraints.
            % Finding the maximum prevalence requires solving a minimization
            % problem, since the resulting equation system involves
            % non-polynomial terms and can not be solved analytically.
            
            highest_degree_monomial = obj.symbolic_penetrances(end);
            bfunction = matlabFunction(rhs(isolate(highest_degree_monomial == 1, sym('b'))));
            lb = 0;
            ub = 1;
            x0 = 0.01;
            nonlcon = @(x) toxo.Model.hconstraint(toxo.PT(obj, maf, x, bfunction(x)), h);
            options = optimoptions('fmincon', 'Display', 'none', 'StepTolerance', 1e-120, 'ConstraintTolerance', 1e-3);
            [alpha, beta, exitflag, output] = fmincon(bfunction, x0, [], [], [], [], lb, ub, nonlcon, options);
            if exitflag < 0
                ME = MException("Model:no_solution", output.message);
                throw(ME);
            end
            pt = toxo.PT(obj, maf, alpha, beta);
        end
        
        function pt = find_max_heritability(obj, maf, p)
            % FIND_MAX_HERITABILITY Compute the penetrance table of the model with the maximum admissible heritability given its MAF and prevalence.
            %   PT = FIND_MAX_PREVALENCE(MAF, H) returns a PT instance with the
            %   maximum admissible heritability given the MAF and prevalence
            %   restraints.
            %   Finding the maximum heritability requires solving a
            %   minimization problem, since the resulting equation system
            %   involves non-polynomial terms and can not be solved
            %   analytically.
            
            highest_degree_monomial = obj.symbolic_penetrances(end);
            afunction = matlabFunction(rhs(isolate(highest_degree_monomial == 1, sym('a'))));
            lb = 0;
            ub = inf;
            x0 = 0;
            nonlcon = @(x) toxo.Model.pconstraint(toxo.PT(obj, maf, afunction(x), x), p);
            options = optimoptions('fmincon', 'Display', 'none', 'StepTolerance', 1e-120, 'ConstraintTolerance', 1e-3);
            [beta, alpha, exitflag, output] = fmincon(afunction, x0, [], [], [], [], lb, ub, nonlcon, options);
            if exitflag < 0
                ME = MException("Model:no_solution", output.message);
                throw(ME);
            end
            pt = toxo.PT(obj, maf, alpha, beta);
        end
    end
end

