function Gamma_out = module5_hermitian_symmetrize(Gamma_in, mode)
% MODULE5_HERMITIAN_SYMMETRIZE - Enforce Hermitian symmetry with real diagonal
%
% Syntax:
%   Gamma_out = module5_hermitian_symmetrize(Gamma_in, mode)
%
% Description:
%   Enforces Hermitian symmetry Γ_ji = conj(Γ_ij) and ensures real diagonal
%   elements. Handles both simplified mode (flexible diagonal) and joint mode
%   (diagonal fixed to 1).
%   
%   Processing steps:
%   1. Compute Hermitian part: (Γ + Γ^H)/2
%   2. Process diagonal according to mode
%   3. Verify final symmetry properties
%
% Input Arguments:
%   Gamma_in - (double array, pxp) Input matrix (potentially non-Hermitian)
%   mode     - (string) Processing mode:
%             'simplified' - Diagonal: real(diagonal)
%             'joint'      - Diagonal: fixed to 1
%
% Output Arguments:
%   Gamma_out - (double array, pxp) Hermitian matrix with real diagonal
%
% Examples:
%   % Basic Hermitian symmetrization
%   Gamma_sym = module5_hermitian_symmetrize(Gamma_rough, 'simplified');
%   
%   % Verify Hermitian property
%   hermitian_error = norm(Gamma_sym - Gamma_sym', 'fro');
%   assert(hermitian_error < 1e-14);
%   
%   % Verify real diagonal
%   diag_imag = max(abs(imag(diag(Gamma_sym))));
%   assert(diag_imag < 1e-14);
%   
%   % Joint mode: diagonal should be exactly 1
%   if strcmp(mode, 'joint')
%       diag_ones_error = norm(diag(Gamma_sym) - ones(size(Gamma_sym, 1), 1));
%       assert(diag_ones_error < 1e-14);
%   end
%
% Mathematical Background:
%   For Hermitian matrices, the symmetry constraint is Γ_ji = conj(Γ_ij).
%   The diagonal must be real for any Hermitian matrix.
%   
%   The symmetrization operation: Γ_sym = (Γ + Γ^H)/2 is the orthogonal
%   projection onto the space of Hermitian matrices in the Frobenius norm.
%
% See also: MODULE5_SINGLE_PROXIMAL_STEP, MODULE5_SOFT_THRESHOLD_COMPLEX
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 2
    error('module5_hermitian_symmetrize:insufficient_input', ...
          'Both Gamma_in and mode are required');
end

% Validate input matrix
if ~isnumeric(Gamma_in) || ~ismatrix(Gamma_in)
    error('module5_hermitian_symmetrize:invalid_input_type', ...
          'Gamma_in must be a numeric matrix');
end

[p, q] = size(Gamma_in);
if p ~= q
    error('module5_hermitian_symmetrize:not_square', ...
          'Gamma_in must be square, got %dx%d', p, q);
end

if ~all(isfinite(Gamma_in(:)))
    error('module5_hermitian_symmetrize:invalid_values', ...
          'Gamma_in must contain finite values');
end

% Validate mode
if ~ischar(mode) && ~isstring(mode)
    error('module5_hermitian_symmetrize:invalid_mode_type', ...
          'mode must be a string');
end

mode = char(mode);
valid_modes = {'simplified', 'joint'};
if ~any(strcmp(mode, valid_modes))
    error('module5_hermitian_symmetrize:invalid_mode_value', ...
          'mode must be one of: %s', strjoin(valid_modes, ', '));
end

% ==================== Hermitian Symmetrization ====================
% Step 1: Project onto Hermitian space via (Γ + Γ^H)/2
Gamma_out = (Gamma_in + Gamma_in') / 2;

% Step 2: Process diagonal according to mode
switch mode
    case 'simplified'
        % Diagonal: take real part, preserve magnitude
        for i = 1:p
            Gamma_out(i, i) = real(Gamma_out(i, i));
        end
        
    case 'joint'
        % Diagonal: fixed to 1 (for joint mode optimization)
        for i = 1:p
            Gamma_out(i, i) = 1.0;
        end
        
    otherwise
        error('module5_hermitian_symmetrize:unsupported_mode', ...
              'Unsupported mode: %s', mode);
end

% ==================== Quality Verification ====================
% Check Hermitian property (should be exact after symmetrization)
hermitian_error = norm(Gamma_out - Gamma_out', 'fro');
if hermitian_error > 1e-12
    warning('module5_hermitian_symmetrize:numerical_error', ...
            'Large Hermitian violation after symmetrization: %.2e', hermitian_error);
end

% Check diagonal reality
diagonal_elements = diag(Gamma_out);
max_diagonal_imag = max(abs(imag(diagonal_elements)));
if max_diagonal_imag > 1e-12
    warning('module5_hermitian_symmetrize:complex_diagonal', ...
            'Large imaginary components in diagonal: %.2e', max_diagonal_imag);
end

% Verify mode-specific diagonal constraints
if strcmp(mode, 'joint')
    diagonal_ones_error = norm(diagonal_elements - ones(p, 1));
    if diagonal_ones_error > 1e-12
        warning('module5_hermitian_symmetrize:diagonal_not_ones', ...
                'Diagonal not exactly ones in joint mode: %.2e', diagonal_ones_error);
    end
end

end