classdef PTable
    %PTable Numeric representation of a penetrance table. This class pro-
    % vides methods to calculate its associated penetrance and heritabi-
    % lity, as well as a method to write the table to a file in several
    % formats.
    
    properties
        order  % Number of loci involved in the penetrance table.
        vars   % Values of the variables present in the original model.
        pt     % Penetrances table array.
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
            %PT Create a penetrance table from a given Model, using the MAF and
            % variable values desired.
            %
            % P = PTable(MODEL, MAFS, VALUES) creates a penetrance table P fol-
            % lowing model description MODEL and using MAFS and variable values
            % inside array VALUES.
            
            obj.order = model.order;
            obj.vars = struct(char(model.variables(1)), values(1), char(model.variables(2)), values(2));
            obj.pt = subs(model.penetrances, model.variables, values);
        end
        
        function p = prevalence(obj, mafs, gp)
            %PREVALENCE Compute the prevalence of the penetrance table.
            if nargin < 3
                gp = toxo.genotype_probabilities(mafs);
            end
            
            p = vpa(sum(obj.pt .* gp));
        end
        
        function h = heritability(obj, mafs)
            %HERITABILITY Compute the heritability of the penetrance table.
            
            gp = toxo.genotype_probabilities(mafs);
            p = obj.prevalence(mafs, gp);
            h = vpa(sum((obj.pt - p).^2 .* gp) / (p * (1 - p)));
        end
        
        function mp = marginal_penetrances(obj, mafs)
            mp = sym(zeros(obj.order, 3));
            rel_penetrances = toxo.genotype_probabilities(mafs) .* obj.pt;
            for p = 1:length(obj.pt)
                for i = 1:obj.order
                    gt = mod(floor((p - 1) / 3^(obj.order - i)), 3) + 1; % index / 3^(order - snp_index) mod 3, with index starting at 0
                    mp(i, gt) = mp(i, gt) + rel_penetrances(p);
                end
            end
            mp = vpa(mp);
        end
        
        function write(obj, path, format, varargin)
            %WRITE Write the penetrance table into a text file using a specific
            % output format.
            %
            % P.WRITE(PATH, FORMAT) writes the penetrance table P into the
            % file specified in PATH using format FORMAT. FORMAT can take any
            % of these values:
            %   -PTable.format_csv: csv-like output file, with each row corre-
            %   sponding to a genotype from the penetrance table and its asso-
            %   ciated phenotype probability.
            %   -PTable.format_gametes: GAMETES compatible model output format.
            
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
