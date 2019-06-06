classdef PTable
    %PTable Symbolic representation of a penetrance table.
    
    properties
        order  % Number of locus defined by the penetrance table.
        vars   % Values for the variables present in the original model.
        pt     % Array of symbolic penetrances values.
    end
    
    methods (Access = private)
        function [s] = to_gametes(obj, fmask, mafs)
            attributes = ['Attribute names:\t' char(join(arrayfun(@(x) sprintf("P%i", x), 0:obj.order-1), char(9))) '\n'];
            solution = char(join(cellfun(@(x) sprintf(string(['%s: ' fmask '\n']), x, vpa(obj.vars.(x))), fieldnames(obj.vars)), ''));
            maf = ['Minor allele frequencies:\t' char(join(arrayfun(@(x) sprintf("%.3f", x), mafs), char(9))) '\n'];
            
            prevalence = sprintf(['Prevalence: ' fmask '\n'], obj.prevalence(mafs));
            heritability = sprintf(['Heritability: ' fmask '\n\n'], obj.heritability(mafs));
            table = ['Table:\n\n' recursive_table(obj.pt, obj.order)];
            s = [attributes maf solution prevalence heritability table];
            
            function s = recursive_table(pt, o)
                n = length(pt) / 3;
                if o > 2
                    s = [recursive_table(pt(1:n), o - 1) '\n' recursive_table(pt(n + 1:2 * n), o - 1) '\n' recursive_table(pt(2 * n + 1:end), o - 1)];
                elseif o == 2
                    s = [recursive_table(pt(1:n), o - 1) recursive_table(pt(n + 1:2 * n), o - 1) recursive_table(pt(2 * n + 1:end), o - 1)];
                else
                    s = [char(join(arrayfun(@(x) string(sprintf(fmask, vpa(x))), pt), ", ")) '\n'];
                end
            end
        end
    end
    
    % Since MATLAB doesn't allow the definition of static variables, they
    % have to be encapsulated inside static methods which can be accessed
    % with no parentheses nor arguments
    methods (Static = true)
        function out = format_csv()
            out = 0;
        end
        
        function out = format_gametes()
            out = 1;
        end
    end
    
    methods
        function obj = PTable(model, values)
            %PTable Create a penetrance table from a given Model and its variable values.
            %
            % PT = PTable(MODEL, VALUES) 
            %   MODEL: Model    Model from which to create the penetrance table.
            %   VALUES: sym     Value for each of the variables represented in MODEL.
            
            obj.order = model.order;
            obj.vars = struct(char(model.variables(1)), values(1), char(model.variables(2)), values(2));
            obj.pt = subs(model.penetrances, model.variables, values);
        end
        
        function p = prevalence(obj, mafs, gp)
            %PREVALENCE Compute the prevalence of the penetrance table.
            
            % P = pt.prevalence(MAFS)
            %   MAFS: double    MAF of each locus.
            %   P: double       Prevalence value.
            
            if nargin < 3
                gp = toxo.genotype_probabilities(mafs);
            end
            
            p = vpa(sum(obj.pt .* gp));
        end
        
        function h = heritability(obj, mafs)
            %HERITABILITY Compute the heritability of the penetrance table.
            %
            % H = pt.heritability(MAFS)
            %   MAFS: double    MAF of each locus.
            %   H: double       Heritability value.
            
            gp = toxo.genotype_probabilities(mafs);
            p = obj.prevalence(mafs, gp);
            h = vpa(sum((obj.pt - p).^2 .* gp) / (p * (1 - p)));
        end
        
        function mp = marginal_penetrances(obj, mafs)
            %MARGINAL_PENETRANCES Compute the marginal penetrance of the three alleles for every locus of the table.
            %
            % MP = pt.marginal_penetrances(MAFS)
            %   MAFS: double    MAF of each locus.
            %   MP: double      Array of the three marginal penetrances for every locus of the table.
            mp = sym(zeros(obj.order, 3));
            gp = toxo.genotype_probabilities(mafs);
            sp = arrayfun(@toxo.genotype_probabilities, mafs, 'UniformOutput', false);
            
            for i = 1:obj.order
                for j = 1:length(obj.pt)
                    geno = mod(fix((j - 1) / 3^(obj.order - i)), 3) + 1;
                    mp(i, geno) = mp(i, geno) + obj.pt(j) * gp(j) / sp{i}(geno);
                end
            end
            mp = vpa(mp);
        end
        
        function write(obj, path, format, varargin)
            %WRITE Write the penetrance table into a file.
            %
            % pt.write(PATH, FORMAT, VARARGIN)
            %   PATH: char      File in which to write the penetrance table.
            %   FORMAT: double  Format to use.
            %   VARARGIN        Additional parameters dependant of the format selected.
            %
            % format_csv: does not require any additional parameters.
            % format_gametes: requires the MAFs of the locus.
            
            fmask = sprintf('%%.%if', fix(digits()/4));
            switch format
                case obj.format_csv
                    s = arrayfun(@(x) {[char(x) char(x)], [char(x) lower(char(x))], [lower(char(x)) lower(char(x))]}, (0:obj.order - 1) + 'A', 'UniformOutput', false);
                    writetable(table(cell2mat(toxo.nfold(s)), arrayfun(@(x) sprintf(fmask, vpa(x)), obj.pt, 'UniformOutput', false)), path, 'WriteVariableNames', false);
                case obj.format_gametes
                    % Variable arguments for GAMETES format:
                    %   mafs: array of doubles
                    fid = fopen(path, 'w+');
                    fprintf(fid, obj.to_gametes(fmask, varargin{1}));
                    fclose(fid);
            end
        end
    end
end
