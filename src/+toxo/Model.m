classdef Model
    %MODEL Symbolic representation of a penetrance model as a function of
    % two variables. This class provides methods to obtain penetrance
    % tables derived from the parametric representation given as input.
    
    properties
        name          % Name of the model.
        order         % Number of loci involved in the epistatic model.
        penetrances   % Array of symbolic expressions, representing the epistatic model.
        variables     % List of all variables contained in all symbolic expressions
    end
    
    methods (Access = private)
        function p = prevalence(obj, maf, gp)
            if nargin < 3
                gp = toxo.genotype_probabilities(maf, obj.order);
            end
            p = sum(obj.penetrances .* gp);
        end
        
        function h = heritability(obj, maf)
            gp = toxo.genotype_probabilities(maf, obj.order);
            p = obj.prevalence(maf, gp);
            h = sum((obj.penetrances - p).^2 .* gp) / (p * (1 - p));
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
            %FIND_MAX_PREVALENCE Calculate the penetrance table(s) of the
            % model with the maximum admissible prevalence given its MAF
            % and heritability.
            %
            % PT = FIND_MAX_PREVALENCE(MAF, H) returns a list of penetrance
            % tables (as PTable objects) that maximize the prevalence, 
            % given the MAF and heritability constraints.
            
            h_constraint = obj.heritability(maf) == h;
            [~, i] = max(subs(obj.penetrances, obj.variables, randi(2^53-1, size(obj.variables))));
            max_poly = obj.penetrances(i);
            [Sx, Sy, p, c] = solve([h_constraint, max_poly == 1], obj.variables, 'Real', true, 'ReturnConditions', true);
            if isempty(Sx)
                ME = MException("Model:no_solution", output.message);
                throw(ME);
            elseif isempty(p)
                pt = toxo.PTable(obj, maf, [Sx, Sy]);
            else
                assume(c);
                Sp = solve(Sx >= 0, Sy >=0, c);
                x = subs(Sx, p, Sp);
                y = subs(Sy, p, Sp);
                pt = toxo.PTable(obj, maf, [x, y]);
            end
        end
        
        function pt = find_max_heritability(obj, maf, p)
            %FIND_MAX_HERITABILITY Calculate the penetrance table(s) of the
            % model with the maximum admissible heritability given its MAF
            % and prevalence.
            % 
            % PT = FIND_MAX_HERITABILITY(MAF, P) returns a list of pene-
            % trance tables (as PTable objects) that maximize the heri-
            % tability, given the MAF and prevalence constraints.
            
            p_constraint = obj.prevalence(maf) == p;
            [~, i] = max(subs(obj.penetrances, obj.variables, randi(2^53-1, size(obj.variables))));
            max_poly = obj.penetrances(i);
            [Sx, Sy, p, c] = solve([p_constraint, max_poly == 1], obj.variables, 'Real', true, 'ReturnConditions', true);
            if isempty(Sx)
                ME = MException("Model:no_solution", output.message);
                throw(ME);
            elseif isempty(p)
                pt = toxo.PTable(obj, maf, [Sx, Sy]);
            else
                assume(c);
                Sp = solve(Sx >= 0, Sy >=0, c);
                x = subs(Sx, p, Sp);
                y = subs(Sy, p, Sp);
                pt = toxo.PTable(obj, maf, [x, y]);
            end
        end
    end
end
