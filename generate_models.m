%% Import Toxo library
addpath('src/');

%% Penetrance table generation
% Read all the models from models/ folder
models = {};
for m = dir('models')'
    if ~ m.isdir
        models{end + 1} = toxo.Model(fullfile(m.folder, m.name));
    end
end

% Create a list of MAFs and heritabilities to test
maf = [0.1, 0.25, 0.4];
h = [0.1, 0.25, 0.5, 0.8];

% Create the output folder
output_folder = "output\";
if ~ isfolder(output_folder)
    mkdir(output_folder);
end

% Find the associated penetrance tables and write the results into files
for m = models
    for i = maf
        for j = h
            try
                pt = m{:}.find_max_prevalence(i, j);
                file_name = sprintf("%s_%.2f_h%.2f.txt", m{:}.name, i, j);
                pt.write(fullfile(output_folder, file_name), toxo.PTable.format_gametes);
            catch ME
                disp(ME.message);
                warning("Unable to generate model %s with MAF %f and heritability %f.\n", m{:}.name, i, j);
                continue;
            end
        end
    end
end
