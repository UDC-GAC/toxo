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
        function m = max_penetrance(obj)
            p = transpose(unique(obj.penetrances(polynomialDegree(obj.penetrances) > 0)));
            m = p(1);
            for e = p(2:end)
                [Sx, ~] = solve([obj.variables >= 0, p >= 0, p <= 1, m > e], 'Real', true);
                if isempty(Sx)
                    m = e;
                end
            end
        end
        
        function [x,y] = solve(obj, varargin)
            assume(in(obj.variables(1), 'Real') & obj.variables(1) >= 0);
            assume(in(obj.variables(2), 'Real') & obj.variables(2) >= 0);
            [Sx, Sy, par, cond] = solve([varargin{:}], obj.variables, 'ReturnConditions', true);
            if isempty(Sx)
                ME = MException("toxo:Model:incompatible", "There is no solution to the problem defined.");
                throwAsCaller(ME);
            elseif isempty(par)
                x = Sx;
                y = Sy;
            else
                % MATLAB does not offer methods to catch warnings, hence it
                % has to be done like this
                lastwarn('');
                warning('off', 'symbolic:solve:SolutionsDependOnConditions');
                warning('off', 'symbolic:solve:ParametrizedFamily');
                S = solve(cond, par);
                warning('on', 'symbolic:solve:SolutionsDependOnConditions');
                warning('on', 'symbolic:solve:ParametrizedFamily');
                [~, wid] = lastwarn();
                if strcmp(wid, 'symbolic:solve:SolutionsDependOnConditions') || strcmp(wid, 'symbolic:solve:ParametrizedFamily')
                    ME = MException("toxo:Model:incapable", "Could not find a solution to the problem defined.");
                    throwAsCaller(ME);
                else
                    x = subs(Sx, par, S);
                    y = subs(Sy, par, S);
                end
            end
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
        
        function pt = find_max_prevalence(obj, mafs, h)
            %FIND_MAX_PREVALENCE Calculate the penetrance table(s) of the
            % model with the maximum admissible prevalence given its MAF
            % and heritability.
            %
            % PT = FIND_MAX_PREVALENCE(MAF, H) returns a list of penetrance
            % tables (as PTable objects) that maximize the prevalence,
            % given the MAF and heritability constraints.
            
            if length(obj.variables) ~= 2
                ME = MException("toxo:Model:unsupported", "Toxo does not support models with %i variables.", length(obj.variables));
                throwAsCaller(ME);
            end
            
            c1 = heritability(mafs) == h;
            c2 = obj.max_penetrance() == 1;
            [x, y] = obj.solve(c1, c2);
            pt = toxo.PTable(obj, [x, y]);
            
            function h = heritability(mafs)
                gp = toxo.genotype_probabilities(mafs);
                p = sum(obj.penetrances .* gp);
                h = sum((obj.penetrances - p).^2 .* gp) / (p * (1 - p));
            end
        end
        
        function pt = find_max_heritability(obj, mafs, p)
            %FIND_MAX_HERITABILITY Calculate the penetrance table(s) of the
            % model with the maximum admissible heritability given its MAF
            % and prevalence.
            %
            % PT = FIND_MAX_HERITABILITY(MAF, P) returns a list of pene-
            % trance tables (as PTable objects) that maximize the heri-
            % tability, given the MAF and prevalence constraints.
            
            if length(obj.variables) ~= 2
                ME = MException("toxo:Model:unsupported", "Toxo does not support models with %i variables.", length(obj.variables));
                throwAsCaller(ME);
            end
            
            c1 = prevalence(mafs) == p;
            c2 = obj.max_penetrance() == 1;
            [x, y] = obj.solve(c1, c2);
            pt = toxo.PTable(obj, [x, y]);
            
            function p = prevalence(mafs)
                gp = toxo.genotype_probabilities(mafs);
                p = sum(obj.penetrances .* gp);
            end
        end
    end
end
