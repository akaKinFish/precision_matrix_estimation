function Gamma_out = module5_active_set_projection(Gamma_in, active_mask)
% MODULE5_ACTIVE_SET_PROJECTION - Project matrix onto active set with Hermitian pairing
%
% Syntax:
%   Gamma_out = module5_active_set_projection(Gamma_in, active_mask)
%
% Description:
%   Projects a matrix onto the active set by zeroing out inactive elements.
%   Ensures Hermitian pairing: if (i,j) is inactive, then (j,i) is also
%   set to zero to maintain symmetry structure.
%   
%   Active set constraint: Γ_ij = 0 for (i,j) ∉ A
%   Hermitian pairing: A(i,j) = A(j,i) always
%
% Input Arguments:
%   Gamma_in    - (double array, pxp) Input matrix
%   active_mask - (logical array, pxp) Active set mask A
%                 true = active element, false = inactive (zero out)
%
% Output Arguments:
%   Gamma_out - (double array, pxp) Matrix with inactive elements zeroed
%
% Examples:
%   % Basic active set projection
%   Gamma_proj = module5_active_set_projection(Gamma, active_mask);
%   
%   % Verify projection properties
%   inactive_elements = find(~active_mask);
%   assert(all(abs(Gamma_proj(inactive_elements)) < 1e-15));
%   
%   % Check Hermitian pairing of active set
%   for i = 1:size(active_mask, 1)
%       for j = 1:size(active_mask, 2)
%           assert(active_mask(i,j) == active_mask(j,i));
%       end
%   end
%   
%   % Verify diagonal is preserved (usually active)
%   if all(diag(active_mask))
%       diag_error = norm(diag(Gamma_proj) - diag(Gamma), 'fro');
%       assert(diag_error < 1e-15);
%   end
%
% Mathematical Background:
%   The active set projection is the orthogonal projection onto the subspace
%   defined by the active set constraints:
%   
%   P_A(Γ)_ij = { Γ_ij  if (i,j) ∈ A
%               { 0     if (i,j) ∉ A
%   
%   For optimization correctness, the active set must respect Hermitian
%   symmetry: (i,j) ∈ A ⟺ (j,i) ∈ A
%
% See also: MODULE5_SINGLE_PROXIMAL_STEP, MODULE3_ACTIVE_SET_MAIN
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 2
    error('module5_active_set_projection:insufficient_input', ...
          'Both Gamma_in and active_mask are required');
end

% Validate input matrix
if ~isnumeric(Gamma_in) || ~ismatrix(Gamma_in)
    error('module5_active_set_projection:invalid_matrix_type', ...
          'Gamma_in must be a numeric matrix');
end

[p, q] = size(Gamma_in);
if p ~= q
    error('module5_active_set_projection:not_square', ...
          'Gamma_in must be square, got %dx%d', p, q);
end

% Validate active mask
if ~islogical(active_mask) || ~ismatrix(active_mask)
    error('module5_active_set_projection:invalid_mask_type', ...
          'active_mask must be a logical matrix');
end

if ~isequal(size(active_mask), [p, p])
    error('module5_active_set_projection:mask_dimension_mismatch', ...
          'active_mask must be %dx%d, got %dx%d', p, p, ...
          size(active_mask, 1), size(active_mask, 2));
end

% ==================== Validate Active Set Properties ====================
% Check Hermitian symmetry of active mask
mask_symmetry_error = sum(sum(active_mask ~= active_mask'));
if mask_symmetry_error > 0
    error('module5_active_set_projection:asymmetric_mask', ...
          'active_mask must be symmetric: A(i,j) = A(j,i). Found %d violations', ...
          mask_symmetry_error);
end

% Check that diagonal is typically active (warning if not)
diagonal_inactive_count = sum(~diag(active_mask));
if diagonal_inactive_count > 0
    warning('module5_active_set_projection:inactive_diagonal', ...
            '%d diagonal elements are inactive, this may cause issues', ...
            diagonal_inactive_count);
end

% ==================== Active Set Projection ====================
% Simple elementwise projection: multiply by mask
Gamma_out = Gamma_in .* active_mask;

% ==================== Verify Projection Properties ====================
% Count active/inactive elements
total_elements = p * p;
active_elements = sum(active_mask(:));
inactive_elements = total_elements - active_elements;

% Verify that inactive elements are indeed zero
inactive_indices = ~active_mask;
max_inactive_value = max(abs(Gamma_out(inactive_indices)));

if max_inactive_value > 1e-15
    warning('module5_active_set_projection:projection_error', ...
            'Some inactive elements not properly zeroed (max: %.2e)', max_inactive_value);
    % Force exact zeros
    Gamma_out(inactive_indices) = 0;
end

% Verify Hermitian symmetry is preserved if input was Hermitian
input_hermitian_error = norm(Gamma_in - Gamma_in', 'fro');
output_hermitian_error = norm(Gamma_out - Gamma_out', 'fro');

if input_hermitian_error < 1e-12 && output_hermitian_error > 1e-12
    warning('module5_active_set_projection:hermitian_violation', ...
            'Projection destroyed Hermitian symmetry (%.2e -> %.2e)', ...
            input_hermitian_error, output_hermitian_error);
end

% ==================== Optional: Store Projection Statistics ====================
% For debugging and monitoring purposes
projection_stats = struct();
projection_stats.total_elements = total_elements;
projection_stats.active_elements = active_elements;
projection_stats.inactive_elements = inactive_elements;
projection_stats.sparsity_ratio = inactive_elements / total_elements;
projection_stats.max_inactive_residual = max_inactive_value;
projection_stats.hermitian_preservation = (output_hermitian_error <= input_hermitian_error + 1e-15);

% Attach stats as persistent property (optional, for advanced debugging)
% This can be retrieved by calling: getappdata(0, 'last_projection_stats')
setappdata(0, 'last_projection_stats', projection_stats);

end