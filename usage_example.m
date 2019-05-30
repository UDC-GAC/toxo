%% Complete usage example
% This script uses the Toxo library to calculate a penetrance table using
% the 2-way additive model using MAF 0.25 and h² 0.2, save it using GAMETES
% format and call GAMETES with that table as input to generate samples.
%
% Author: Christian Ponte
% Last modified: 30 May 2019

%% Script
% Import Toxo library
addpath('src/');

% Path definitions
gametes = '~/Tools/Epistasis/Simulation/GAMETES/GAMETES_2.1.jar';

% Model definitions
model = 'models/additive_2.csv';
mafs = [0.25, 0.25];
h = 0.2;

% Create the output folder
output_folder = "output";
if ~ isfolder(output_folder)
    mkdir(output_folder);
end

%%% TRY

% Create a Model instance
m = toxo.Model(model);
% Calculate the penetrance table with maximum P(D)
pt = m.find_max_heritability(mafs, h);
% Save the table using GAMETES format
table_path = char(fullfile(output_folder, 'table.txt'));
pt.write(table_path, toxo.PTable.format_gametes, mafs);

% Call GAMETES to generate samples
data_path = char(fullfile(output_folder, 'example'));
system(['java -jar ' gametes ' -i ' table_path ' -D "-n 0.05 -x 0.5 -a 1000 -s 100 -w 100 -r 10 -o ' data_path '"']);
