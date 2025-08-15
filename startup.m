%% Startup script for Precision Matrix Estimation project
% This script sets up the MATLAB path for the project

% Get the directory of this script
project_root = fileparts(mfilename('fullpath'));
if isempty(project_root)
    project_root = pwd;
end

% Add all necessary directories to path
addpath(project_root);
addpath(genpath(fullfile(project_root, 'src')));
addpath(genpath(fullfile(project_root, 'utils')));
addpath(genpath(fullfile(project_root, 'tests')));
addpath(fullfile(project_root, 'examples'));
addpath(fullfile(project_root, 'config'));

fprintf('Precision Matrix Estimation project paths configured.\\n');
fprintf('Project root: %s\\n', project_root);
