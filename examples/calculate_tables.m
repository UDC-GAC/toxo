%% Table generation script
% This script uses the Toxo library to calculate penetrance tables of the
% additive, multiplicative and threshold models, for third and fourth
% order, using MAFs 0.1 and 0.4 for all locus and heritabilities 0.1 and
% 0.8.
%
% Author: Christian Ponte
% Last modified: 30 May 2019

%% Script
% Import Toxo library
addpath('src/');

% Read all the third and fourth order models from models/ folder
model_list = {
    'models/additive_3.csv'
    'models/multiplicative_3.csv'
    'models/threshold_3.csv'
    'models/additive_4.csv'
    'multiplicative_4.csv'
    'threshold_4.csv'
    };

% Create a list of MAFs and heritabilities to test
maf = [0.1, 0.4];
h = [0.1, 0.8];

% Create the output folder
output_folder = "output";
if ~ isfolder(output_folder)
    mkdir(output_folder);
end

for i = 1:length(model_list)
    % Create a Model instance
    m = toxo.Model(model_list{i});
    for j = 1:length(maf)
        mafs = repmat(maf(j), 1, m.order);
        for k = 1:length(h)
            try
                % Calculate a penetrance table for the model
                pt = m.find_max_prevalence(mafs, h(k));
                % Write the table to a file
                file_name = sprintf("%s_%.1f_h%.1f.txt", m.name, maf(j), h(k));
                pt.write(fullfile(output_folder, file_name), toxo.PTable.format_csv);
            catch ME
                warning("(%s with MAF=%.1f and h²=%.1f) %s", m.name, maf(j), h(k), ME.message);
                continue;
            end
        end
    end
end
