classdef PTable
    %PTable Numeric representation of a penetrance table. This class pro-
    % vides methods to calculate its associated penetrance and heritabi-
    % lity, as well as a method to write the table to a file in several
    % formats.
    
    properties
        order  % Number of loci involved in the penetrance table.
        maf    % Common MAF of all locis involved in the interaction.
        val    % Values of the two variables in the original model.
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
                s = join(arrayfun(@(x) sprintf("%.12f", x), pt), ", ") + newline;
            end
        end
    end
    
    % Since MATLAB doesn't allow the definition of static variables, they have to be encapsulated inside static methods with the "persistent" definition
    methods (Static = true)
        function out = format_plaintext()
            persistent plaintext;
            if isempty(plaintext)
                plaintext = toxo.PTable.counter;
            end
            out = plaintext;
        end
        
        function out = format_gametes()
            persistent gametes;
            if isempty(gametes)
                gametes = toxo.PTable.counter;
            end
            out = gametes;
        end
    end
    
    methods
        function obj = PTable(model, maf, values)
        % PT Create a penetrance table from a given Model, using the MAF, alpha and beta desired.
        %   P = PT(MODEL, M, A, B) creates a penetrance table P following
        %   model description MODEL and using MAF M, baseline effect A and
        %   genotypic effect B.
            
            obj.order = model.order;
            obj.maf = maf;
            obj.val = values;
            obj.gp = toxo.genotype_probabilities(maf, model.order);
            obj.pt = double(subs(model.symbolic_penetrances, model.variables, [x, y]));
        end
        
        function p = prevalence(obj)
        % PREVALENCE Compute the prevalence of the penetrance table.

            p = sum(obj.pt .* obj.gt_p);
        end
        
        function h = heritability(obj)
        % HERITABILITY Compute the heritability of the penetrance table.

            p = obj.prevalence();
            h = sum((obj.pt - p).^2 .* obj.gt_p) / p * (1 - p);
        end
        
        function write(obj, path, format)
        % WRITE Write the penetrance into a text file using a specific output format.
        %   P.WRITE(PATH, FORMAT) writes the penetrance table P into the
        %   file specified in PATH using format FORMAT. FORMAT can take any
        %   of these values:
        %     -PT.format_paintext: text-format output file, with each row
        %     corresponding to a genotype from the penetrance table and its
        %     associated phenotype probability.
        %     -PT.format_gametes: GAMETES compatible model output format.
        
            switch format
                case obj.format_plaintext
                    gt_strings = char(char(join(toxo.nfold(["AA" "Aa" "aa"], obj.order), "")) + repelem(0:obj.order - 1, 2));
                    fid = fopen(path, 'w+');
                    for i = 1:length(obj.pt)
                        fprintf(fid, "%s %.12f\n", gt_strings(i, :), obj.pt(i));
                    end
                    fclose(fid);
                case obj.format_gametes
                    header_template = ...
                        "Attribute names:\t%s\n" + ...
                        "Minor allele frequencies:\t%s\n" + ...
                        "Alpha: %.12f\n" + ...
                        "Beta: %.12f\n" + ...
                        "K: %f\n" + ...
                        "Heritability: %f\n\n" + ...
                        "Table:\n\n" + ...
                        "%s";
                    names = join(arrayfun(@(x) sprintf("P%i", x), 0:obj.order-1), char(9));
                    maf_list = join(arrayfun(@(x) sprintf("%.3f", x), repmat(obj.maf, [1, obj.order])), char(9));
                    p = obj.prevalence();
                    h = obj.heritability();
                    fid = fopen(path, 'w+');
                    fprintf(fid, header_template, names, maf_list, obj.alpha, obj.beta, p, h, toxo.PTable.pt_to_string_table(obj.pt, obj.order));
                    fclose(fid);
            end
        end
    end
end
