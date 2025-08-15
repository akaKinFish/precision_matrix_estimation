function params_struct = convert_params_to_struct(varargin)
% CONVERT_PARAMS_TO_STRUCT - Convert name-value pairs to struct
% Location: utils/convert_params_to_struct.m
%
% This helper function safely converts parameter name-value pairs 
% (either as cell array or individual arguments) to a struct.
%
% Usage:
%   params = convert_params_to_struct('n_nodes', 10, 'n_freq', 20);
%   params = convert_params_to_struct(param_cell{:});
%   params = convert_params_to_struct(param_cell);  % Also handles cell array directly
%
% Example:
%   common_params = {'n_nodes', 15, 'n_freq', 40};
%   params = convert_params_to_struct(common_params);
%   [prec, cov, emp, p] = module7_simulation_from_struct(params);

    if nargin == 1 && iscell(varargin{1})
        % Handle case where a cell array is passed directly
        param_cell = varargin{1};
    else
        % Handle case where parameters are passed as varargin
        param_cell = varargin;
    end
    
    % Ensure we have an even number of elements
    if mod(length(param_cell), 2) ~= 0
        error('Parameters must be in name-value pairs');
    end
    
    % Convert to struct
    params_struct = struct();
    for i = 1:2:length(param_cell)
        param_name = param_cell{i};
        param_value = param_cell{i+1};
        
        % Validate parameter name
        if ~ischar(param_name) && ~isstring(param_name)
            error('Parameter name must be a string at position %d', i);
        end
        
        % Store in struct
        params_struct.(param_name) = param_value;
    end
end