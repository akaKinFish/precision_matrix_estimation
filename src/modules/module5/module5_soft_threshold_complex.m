function Gamma_out = module5_soft_threshold_complex(Gamma_in, tau, mode)
% MODULE5_SOFT_THRESHOLD_COMPLEX - Complex amplitude soft thresholding with Hermitian symmetry
%
% Syntax:
%   Gamma_out = module5_soft_threshold_complex(Gamma_in, tau, mode)
%
% Description:
%   Applies complex amplitude soft thresholding to off-diagonal elements while
%   preserving phase information and ensuring Hermitian symmetry. The proximal
%   operator for L1 penalty on complex matrices.
%   
%   For complex z: prox_τ|·|(z) = (1 - τ/|z|)z if |z| > τ, 0 otherwise
%   
%   Processing order:
%   1. Apply thresholding to upper triangular elements (i < j)
%   2. Mirror to lower triangular with conjugate symmetry
%   3. Process diagonal according to mode
%
% Input Arguments:
%   Gamma_in - (double array, pxp) Input matrix (may be complex)
%   tau      - (double) Thresholding parameter τ = λ₂α
%   mode     - (string) Processing mode:
%             'simplified' - Diagonal: real(diagonal), no thresholding
%             'joint'      - Diagonal: fixed to 1
%
% Output Arguments:
%   Gamma_out - (double array, pxp) Thresholded matrix with Hermitian symmetry
%
% Examples:
%   % Basic soft thresholding
%   Gamma_thresh = module5_soft_threshold_complex(Gamma, 0.01, 'simplified');
%   
%   % Verify Hermitian property
%   hermitian_error = norm(Gamma_thresh - Gamma_thresh', 'fro');
%   assert(hermitian_error < 1e-12);
%   
%   % Check phase preservation for non-zero elements
%   [i, j] = find(abs(Gamma_thresh) > 1e-12);
%   for k = 1:length(i)
%       if i(k) ~= j(k)  % Off-diagonal
%           original_phase = angle(Gamma_in(i(k), j(k)));
%           new_phase = angle(Gamma_thresh(i(k), j(k)));
%           phase_diff = abs(original_phase - new_phase);
%           assert(phase_diff < 1e-12 || phase_diff > 2*pi - 1e-12);
%       end
%   end
%
% Mathematical Background:
%   The soft thresholding operator for complex numbers preserves the phase
%   while shrinking the amplitude:
%   
%   If |z| > τ: result = |z| * (|z| - τ)/|z| * exp(i*angle(z)) = (1 - τ/|z|)*z
%   If |z| ≤ τ: result = 0
%   
%   For Hermitian matrices: Γ_ji = conj(Γ_ij), so we process upper triangle
%   and mirror with conjugate to ensure proper symmetry.
%
% See also: MODULE5_SINGLE_PROXIMAL_STEP, MODULE5_HERMITIAN_SYMMETRIZE
%
% Author: [Your Name]
% Date: [Current Date]
% Version: 1.0

% ==================== Input Validation ====================
if nargin < 3
    error('module5_soft_threshold_complex:insufficient_input', ...
          'All 3 arguments are required');
end

% Validate input matrix
if ~isnumeric(Gamma_in) || ~ismatrix(Gamma_in)
    error('module5_soft_threshold_complex:invalid_input_type', ...
          'Gamma_in must be a numeric matrix');
end

[p, q] = size(Gamma_in);
if p ~= q
    error('module5_soft_threshold_complex:not_square', ...
          'Gamma_in must be square, got %dx%d', p, q);
end

% Validate threshold parameter
if ~isscalar(tau) || ~isreal(tau) || tau < 0
    error('module5_soft_threshold_complex:invalid_threshold', ...
          'tau must be a non-negative real scalar');
end

% Validate mode
if ~ischar(mode) && ~isstring(mode)
    error('module5_soft_threshold_complex:invalid_mode', ...
          'mode must be a string');
end

mode = char(mode);
valid_modes = {'simplified', 'joint'};
if ~any(strcmp(mode, valid_modes))
    error('module5_soft_threshold_complex:invalid_mode_value', ...
          'mode must be one of: %s', strjoin(valid_modes, ', '));
end

% ==================== Complex Soft Thresholding ====================
% Initialize output matrix
Gamma_out = Gamma_in;

% Process upper triangular elements (i < j) with complex amplitude thresholding
for i = 1:p
    for j = i+1:p  % Upper triangle only
        z = Gamma_in(i, j);
        magnitude = abs(z);
        
        if magnitude > tau && magnitude > 1e-15  % Avoid division by very small numbers
            % Shrink amplitude while preserving phase
            shrinkage_factor = 1 - tau / magnitude;
            Gamma_out(i, j) = shrinkage_factor * z;
        else
            % Threshold to zero
            Gamma_out(i, j) = 0;
        end
    end
end

% Mirror to lower triangular with conjugate symmetry: Γ_ji = conj(Γ_ij)
for i = 1:p
    for j = 1:i-1  % Lower triangle
        Gamma_out(i, j) = conj(Gamma_out(j, i));
    end
end

% ==================== Diagonal Processing ====================
switch mode
    case 'simplified'
        % Diagonal: take real part, no thresholding
        for i = 1:p
            Gamma_out(i, i) = real(Gamma_in(i, i));
        end
        
    case 'joint'
        % Diagonal: fixed to 1 (for joint mode with Δ updates)
        for i = 1:p
            Gamma_out(i, i) = 1.0;
        end
        
    otherwise
        error('module5_soft_threshold_complex:unsupported_mode', ...
              'Unsupported mode: %s', mode);
end

% ==================== Final Hermitian Verification ====================
% Ensure numerical Hermitian symmetry
hermitian_error = norm(Gamma_out - Gamma_out', 'fro');
if hermitian_error > 1e-12
    % Force perfect symmetry
    Gamma_out = (Gamma_out + Gamma_out') / 2;
    
    % Re-ensure real diagonal
    for i = 1:p
        if strcmp(mode, 'simplified')
            Gamma_out(i, i) = real(Gamma_out(i, i));
        else  % joint mode
            Gamma_out(i, i) = 1.0;
        end
    end
end

% ==================== Quality Checks ====================
% Verify final properties (for debugging)
final_hermitian_error = norm(Gamma_out - Gamma_out', 'fro');
final_diagonal_imag = max(abs(imag(diag(Gamma_out))));

if final_hermitian_error > 1e-10
    warning('module5_soft_threshold_complex:hermitian_violation', ...
            'Large Hermitian violation after processing: %.2e', final_hermitian_error);
end

if final_diagonal_imag > 1e-10
    warning('module5_soft_threshold_complex:complex_diagonal', ...
            'Large imaginary diagonal components: %.2e', final_diagonal_imag);
end

end