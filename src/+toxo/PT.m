classdef PT
    %PT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        format_plaintext = 0
        format_gametes = 1
    end
    
    properties
        order
        maf
        alpha
        beta
        gt_p
        pt
    end
    
    methods (Static = true, Access = protected)
        function [s] = pt_to_string_table(pt, o)
            n = length(pt) / 3;
            if o > 2
                s = penetrance_list_to_string(pt(1:n), o - 1) + newline + ...
                    penetrance_list_to_string(pt(n + 1:2 * n), o - 1) + ...
                    newline + penetrance_list_to_string(pt(2 * n + 1:end), o - 1);
            elseif o == 2
                s = penetrance_list_to_string(pt(1:n), o - 1) + ...
                    penetrance_list_to_string(pt(n + 1:2 * n), o - 1) + ...
                    penetrance_list_to_string(pt(2 * n + 1:end), o - 1);
            else
                s = join(arrayfun(@(x) sprintf("%f", x), pt), ", ") + newline;
            end
        end
    end
    
    methods
        function obj = PT(model, maf, alpha, beta)
            %PT Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.order = model.order;
            obj.maf = maf;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.gt_p = model.genotype_probabilities(maf);
            obj.pt = double(subs( ...
                model.pt_s, [sym('a'), sym('b')], [alpha, beta] ...
            ));
        end
        
        function p = prevalence(obj)
            % PREVALENCE Compute the prevalence of the given penetrance table using
            %   the genotype combination probabilities.

            p = sum(obj.pt .* obj.gt_p);
        end
        
        function h = heritability(obj)
            % HERITABILITY Compute the heritability of the given penetrance table using
            %   its prevalence value and the genotype combination probabilities.

            p = obj.prevalence();
            h = sum((obj.pt - p).^2 .* obj.gt_p) / p * (1 - p);
        end
        
        function write(obj, path, format)
            switch format
                case obj.format_plaintext
                    gt_strings = toxo.nfold(["AA" "Aa" "aa"], obj.order);
                    fid = fopen(path, 'w+');
                    for i = 1:length(obj.pt)
                        fprintf(fid, "%s %.2f\n", gt_strings(i, :), obj.pt(i));
                    end
                    fclose(fid);
                case obj.format_gametes
                    header_template = ...
                        "Attribute names:\t%s\n" + ...
                        "Minor allele frequencies:\t%s\n" + ...
                        "Alpha: %f\n" + ...
                        "Beta: %f\n" + ...
                        "K: %f\n" + ...
                        "Heritability: %f\n\n" + ...
                        "Table:\n\n" + ...
                        "%s";
                    names = join(arrayfun(@(x) sprintf("P%i", x), ...
                        0:obj.order-1), "	");
                    maf_list = join(arrayfun(@(x) sprintf("%.2f", x), ...
                        repmat(obj.maf, [1, obj.order])), "	");
                    p = obj.prevalence();
                    h = obj.heritability();
                    fid = fopen(path, 'w+');
                    fprintf(fid, header_template, names, maf_list, ...
                        obj.alpha, obj.beta, p, h, ...
                        pt_to_string_table(obj.pt, obj.order));
                    fclose(fid);
            end
        end
    end
end
