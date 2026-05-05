function P = omega_to_pcor(Omega, diag_mode)
    % Omega: 一个频率点的 precision matrix (Nr x Nr), 可为 complex Hermitian
    % diag_mode: 'one' | 'zero' | 'nan'
    if nargin < 2
        diag_mode = 'zero';
    end

    % 数值上先强制 Hermitian
    G = (Omega + Omega') / 2;

    % 分母只用实对角，防止除 0
    d = max(real(diag(G)), 1e-12);

    % partial correlation / partial coherence 风格归一化
    P = -G ./ sqrt(d * d.');

    % 对角线怎么放，按你的用途来
    switch lower(diag_mode)
        case 'one'
            P(1:size(G,1)+1:end) = 1;
        case 'zero'
            P(1:size(G,1)+1:end) = 0;
        case 'nan'
            P(1:size(G,1)+1:end) = NaN;
        otherwise
            error('diag_mode must be one / zero / nan');
    end
end