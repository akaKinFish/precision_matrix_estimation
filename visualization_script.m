opts = struct();
opts.alpha_band = [8 12];
% ROI chord (很多边) → quantile 高一点
opts.roi_edge_keep_quantile = 5;
opts.roi_edge_keep_max      = 150;
% Region chord (很少边) → quantile 必须低很多
opts.region_edge_keep_quantile = 5;
opts.region_edge_keep_max      = 150;
opts.region_stat = 'both';
opts.debug_print = true;
opts.debug_topK  = 10;
opts.SC  = Sc;
opts.atlas_name = 'HCP_MMP1';
opts.save_dir = 'F:\hope_the_final\figs_G5IWCMKRGRXQ_v2';
out = visualize_subject_jspace_result_v2( ...
"D:\results\sub11_jspace_results.mat", ...
"D:\result\Switzerland\2pre\crossSpec\G5IWCMKRGRXQ\G5IWCMKRGRXQ.mat", opts);