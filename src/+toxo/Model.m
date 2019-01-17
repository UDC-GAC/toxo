classdef Model
    %MODEL Abstract representation of a penetrance model as a function of
    % two variables. This class provides methods to obtain numeric pene-
    % trance tables derived from the parametric representation given as
    % input.
    
    properties
        name          % Name of the model.
        order         % Number of loci involved in the epistatic model.
        penetrances   % Array of symbolic expressions, representing the epistatic model.
        variables     % List of all variables contained in all symbolic expressions
    end
    
    methods (Access = private)
        function p = prevalence(obj, maf)
            gp = toxo.genotype_probabilities(maf, obj.order);
            p = sum(obj.penetrances .* gp);
        end
    end
    
    methods
        function obj = Model(path)
            %MODEL Construct an instance of this class from the given
            % model.
            % 
            % M = MODEL(P) reads the model from its text representation in
            % file path P. The input model must be formatted as a plain
            % CSV, with eachline of the file P corresponding to a row of 
            % the model. The rows are made of the genotype definition and
            % the probability associated with the given genotype, separated
            % by a coma. Probability is expressed as a function of two 
            % variables, represented as alphabetical characters. Empty
            % lines, as well as lines starting with '#' will be ignored.
            
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
            % FIND_MAX_PREVALENCE Compute the penetrance table of the model with the maximum admissible prevalence given its MAF and heritability.
            %   PT = FIND_MAX_PREVALENCE(MAF, H) returns a PT instance with the
            %   maximum admissible prevalence given the MAF and heritability
            %   restraints.
            %   Finding the maximum prevalence requires solving a minimization
            %   problem, since the resulting equation system involves
            %   non-polynomial terms and can not be solved analytically.
            
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
            %FIND_MAX_HERITABILITY Calculate the penetrance table(s) of the
            % model with the maximum admissible heritability given its MAF
            % and prevalence.
            % 
            % PT = FIND_MAX_HERITABILITY(MAF, P) returns a list of pene-
            % trance tables (as PTable objects) that maximize the heri-
            % tability given the MAF and prevalence constraints.
            
            p_constraint = obj.prevalence(maf) == p;
            [~, i] = max(subs(obj.penetrances, obj.variables, randi(2^53-1, size(obj.variables))));
            max_poly = obj.penetrances(i);
            S = solve([p_constraint, max_poly == 1], obj.variables);
            x = vpa(S.x);
            y = vpa(S.y);
            solutions = arrayfun(@(z, t) isreal(z) && isreal(t) && z >= 0 && t >= 0, x, y);
            x = x(solutions);
            y = y(solutions);
            pt = arrayfun(@(z, t) toxo.PT(obj, maf, z, t), x, y);
        end
    end
end

