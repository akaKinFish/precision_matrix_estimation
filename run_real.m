data = load('G:\BC-V_result_12_6_part1\CHBMP\2pre\crossSpec\F5GTSSCCIRIQ\F5GTSSCCIRIQ.mat');

Svv_CrossM = data.data_struct.CrossM;
dnames = data.data_struct.dnames;
parameters = load('G:\OneDrive - CCLAB\New_Data_XiAlphaNET\TestRetest_data\preprocessed data\session1\structural\parameters.mat');

C = parameters.Model.C;
L = parameters.Model.K;

compact_C = parameters.Compact_Model.C;
compact_L = parameters.Compact_Model.K;
freq = data.data_struct.freqrange;

[Omega_est, Sjj_est, outs_js] = run_jspace_real(Svv_CrossM, compact_L, freq, compact_C);

% Omega_est, Sjj_est: 47x1 cell, each cell: 360x360 complex double
Nw   = numel(Sjj_est);
Nroi = size(Sjj_est{1},1);

Omega3 = zeros(Nroi, Nroi, Nw, 'like', Omega_est{1});
Sjj3   = zeros(Nroi, Nroi, Nw, 'like', Sjj_est{1});

for k = 1:Nw
    Omega3(:,:,k) = Omega_est{k};
    Sjj3(:,:,k)   = Sjj_est{k};

    % 可选：数值安全（Hermitian + real diagonal）
    A = Omega3(:,:,k); A = (A + A')/2; A(1:Nroi+1:end) = real(diag(A)); Omega3(:,:,k) = A;
    B = Sjj3(:,:,k);   B = (B + B')/2; B(1:Nroi+1:end) = real(diag(B)); Sjj3(:,:,k)   = B;
end

freq = freq(:);                 % 47x1
Nw   = length(freq);
Nroi = size(Sjj3,1);

pow = zeros(Nroi, Nw);
for k = 1:Nw
    dk = real(diag(Sjj3(:,:,k)));
    dk = max(dk, 0);            % 防止数值误差出现负
    pow(:,k) = dk;
end
logPow = log10(pow + eps);      % Nroi x Nw

figure('Color','w'); 
plot(freq, logPow, 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Mean log_{10}(Power)');
title('Mean Source Log Spectrum (diag(Sjj))'); grid on;

roi_id = 1;
figure('Color','w');
plot(freq, logPow(roi_id,:), 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel(sprintf('ROI %d log_{10}(Power)', roi_id));
title(sprintf('ROI %d Log Spectrum', roi_id)); grid on;

q25 = prctile(logPow, 25, 1);
q50 = prctile(logPow, 50, 1);
q75 = prctile(logPow, 75, 1);

figure('Color','w'); hold on;
fill([freq; flipud(freq)], [q25(:); flipud(q75(:))], 0.8*[1 1 1], 'EdgeColor','none');
plot(freq, q50, 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('log_{10}(Power)');
title('Source Log Spectrum: median + IQR across ROIs'); grid on;

Nroi = size(Sjj3,1);
Nw   = size(Sjj3,3);

Coh3 = zeros(Nroi, Nroi, Nw);
for k = 1:Nw
    S = Sjj3(:,:,k);
    p = real(diag(S)); 
    p = max(p, 0) + eps;
    denom = p * p.';                 % (Sii*Sjj)
    Coh3(:,:,k) = (abs(S).^2) ./ denom;
    Coh3(1:Nroi+1:end,k) = 1;        % 可选：对角线强制1（或者你也可以置0）
end

meanCoh = zeros(Nw,1);
mask = ~eye(Nroi);

for k = 1:Nw
    C = Coh3(:,:,k);
    meanCoh(k) = mean(C(mask), 'omitnan');
end

figure('Color','w');
plot(freq, meanCoh, 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Mean Coherence (off-diagonal)');
title('Global Mean Coherence vs Frequency'); grid on;


idx = find(freq>=8 & freq<=13);
C_alpha = mean(Coh3(:,:,idx), 3, 'omitnan');
C_alpha(1:Nroi+1:end) = 0;   % 热图通常把对角线去掉

figure('Color','w');
imagesc(C_alpha); axis image; colorbar;
title('Alpha-band Coherence (mean over 8-13 Hz)');
xlabel('ROI'); ylabel('ROI');

