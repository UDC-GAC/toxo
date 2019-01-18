classdef PTable
    %PTable Numeric representation of a penetrance table. This class pro-
    % vides methods to calculate its associated penetrance and heritabi-
    % lity, as well as a method to write the table to a file in several
    % formats.
    
    properties
        order  % Number of loci involved in the penetrance table.
        maf    % Common MAF of all locis involved in the interaction.
        vars   % Values of the two variables in the original model.
        gp     % Genotype probability table array.
        pt     % Penetrance table array.
    end
    
    methods (Access = private, Static = true)
        function n = counter()
            persistent value;
            if isempty(value)
                value = 0;
            end
            n = value;
            value = value + 1;
        end
        
        function [s] = pt_to_string_table(pt, o)
            n = length(pt) / 3;
            if o > 2
                s = toxo.PTable.pt_to_string_table(pt(1:n), o - 1) + newline + ...
                    toxo.PTable.pt_to_string_table(pt(n + 1:2 * n), o - 1) + newline + ...
                    toxo.PTable.pt_to_string_table(pt(2 * n + 1:end), o - 1);
            elseif o == 2
                s = toxo.PTable.pt_to_string_table(pt(1:n), o - 1) + ...
                    toxo.PTable.pt_to_string_table(pt(n + 1:2 * n), o - 1) + ...
                    toxo.PTable.pt_to_string_table(pt(2 * n + 1:end), o - 1);
            else
                s = join(arrayfun(@(x) sprintf("%.12f", vpa(x)), pt), ", ") + newline;
            end
        end
    end
    
    % Since MATLAB doesn't allow the definition of static variables, they have to be encapsulated inside static methods with the "persistent" definition
    methods (Static = true)
        function out = format_csv()
            out = 0;
        end
        
        function out = format_gametes()
            out = 1;
        end
    end
    
    methods
        function obj = PTable(model, maf, values)
        %PT Create a penetrance table from a given Model, using the MAF and variable values desired.
        % P = PTable(MODEL, M, A, B) creates a penetrance table P following
        % model description MODEL and using MAF M, baseline effect A and
        % genotypic effect B.
            
            obj.order = model.order;
            obj.maf = maf;
            obj.vars = containers.Map(arrayfun(@char, model.variables, 'uniform', 0), double(values));
            obj.gp = toxo.genotype_probabilities(maf, model.order);
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
        %WRITE Write the penetrance table into a text file using a specific output format.
        % P.WRITE(PATH, FORMAT) writes the penetrance table P into the
        % file specified in PATH using format FORMAT. FORMAT can take any
        % of these values:
        %   -PTable.format_paintext: text-format output file, with each row
        %   corresponding to a genotype from the penetrance table and its
        %   associated phenotype probability.
        %   -PTable.format_gametes: GAMETES compatible model output format.
        
            switch format
                case obj.format_csv
                    gt_strings = string(char(char(join(toxo.nfold(["AA" "Aa" "aa"], obj.order), "")) + repelem(0:obj.order - 1, 2)));
                    writetable(table(gt_strings, arrayfun(@(x) sprintf("%.12f", vpa(x)), obj.pt)), path, 'WriteVariableNames', false);
                case obj.format_gametes
                    header_template = ...
                        "Attribute names:\t%s\n" + ...
                        "Minor allele frequencies:\t%s\n" + ...
                        "%s" + ...
                        "Prevalence: %f\n" + ...
                        "Heritability: %f\n\n" + ...
                        "Table:\n\n" + ...
                        "%s";
                    names = join(arrayfun(@(x) sprintf("P%i", x), 0:obj.order-1), char(9));
                    maf_list = join(arrayfun(@(x) sprintf("%.3f", x), repmat(obj.maf, [1, obj.order])), char(9));
                    vars_list = join(cellfun(@(x) sprintf("%s: %.12f\n", x, obj.vars(x)), obj.vars.keys()), "");
                    fid = fopen(path, 'w+');
                    fprintf(fid, header_template, names, maf_list, vars_list, obj.prevalence, obj.heritability, toxo.PTable.pt_to_string_table(obj.pt, obj.order));
                    fclose(fid);
            end
        end
    end
end
