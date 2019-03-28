classdef PTable
    %PTable Numeric representation of a penetrance table. This class pro-
    % vides methods to calculate its associated penetrance and heritabi-
    % lity, as well as a method to write the table to a file in several
    % formats.
    
    properties
        order  % Number of loci involved in the penetrance table.
        mafs   % MAF associated with each loci involved in the interaction.
        vars   % Values of the variables present in the original model.
        gp     % Genotype probabilities table array.
        pt     % Penetrances table array.
    end
    
    methods (Access = private)
        function [s] = to_gametes(obj, fmask)
            attributes = ['Attribute names:\t' char(join(arrayfun(@(x) sprintf("P%i", x), 0:obj.order-1), char(9))) '\n'];
            maf = ['Minor allele frequencies:\t' char(join(arrayfun(@(x) sprintf("%.3f", x), obj.mafs), char(9))) '\n'];
            solution = char(join(cellfun(@(x) sprintf(string(['%s: ' fmask '\n']), x, vpa(obj.vars(x))), obj.vars.keys()), ''));
            prevalence = sprintf(['Prevalence: ' fmask '\n'], obj.prevalence);
            heritability = sprintf(['Heritability: ' fmask '\n\n'], obj.heritability);
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
        function obj = PTable(model, mafs, values)
        %PT Create a penetrance table from a given Model, using the MAF and
        % variable values desired.
        % 
        % P = PTable(MODEL, MAFS, VALUES) creates a penetrance table P fol-
        % lowing model description MODEL and using MAFS and variable values
        % inside array VALUES.
            
            obj.order = model.order;
            obj.mafs = mafs;
            obj.vars = containers.Map(arrayfun(@char, model.variables, 'uniform', 0), arrayfun(@(x) {x}, values));
            obj.gp = toxo.genotype_probabilities(mafs);
            obj.pt = subs(model.penetrances, model.variables, values);
        end
        
        function p = prevalence(obj)
        %PREVALENCE Compute the prevalence of the penetrance table.

            p = vpa(sum(obj.pt .* obj.gp));
        end
        
        function h = heritability(obj)
        %HERITABILITY Compute the heritability of the penetrance table.

            p = sum(obj.pt .* obj.gp);
            h = vpa(sum((obj.pt - p).^2 .* obj.gp) / (p * (1 - p)));
        end
        
        function write(obj, path, format)
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
                    s = arrayfun(@(x) {[char(x) char(x)], [char(x) lower(char(x))], [lower(char(x)) lower(char(x))]}, (0:length(obj.mafs) - 1) + 'A', 'UniformOutput', false);
                    writetable(table(cell2mat(toxo.nfold(s)), arrayfun(@(x) sprintf(fmask, vpa(x)), obj.pt, 'UniformOutput', false)), path, 'WriteVariableNames', false);
                case obj.format_gametes
                    fid = fopen(path, 'w+');
                    fprintf(fid, obj.to_gametes(fmask));
                    fclose(fid);
            end
        end
    end
end
