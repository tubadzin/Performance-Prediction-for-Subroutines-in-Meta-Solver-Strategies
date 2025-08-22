function plot_model_perf_capslack(scen,capping_type)
% load(['fixed_thresh_results/matrix/',scen,'-cens-matrix-pred-scale']);
% load(['results/matrix/',scen,'/a_cens-matrix-pred-scale']); % saved values in there for capslack
% load(['results/matrix/',scen,'/cens-matrix-pred-scale-capslack']);

% capping_type = 'fixed';
switch capping_type
    case 'fixed'
        load(['results/matrix/',scen,'/cens-matrix-pred-scale-fixed']);
        xLabelString = 'Fixed censoring threshold [s]';
    case 'capslack'
        load(['results/matrix/',scen,'/cens-matrix-pred-scale-capslack']);
        xLabelString = 'Slack factor for censoring';
    otherwise
        error
end

for cs_idx=1:length(caps)
    for strategy = 1:4
        for seed=1:1
            rmses(seed) = results{strategy,cs_idx,seed}{1}(1); 
            lls(seed) = results{strategy,cs_idx,seed}{1}(2); 
            ccs(seed) = results{strategy,cs_idx,seed}{1}(3); 
        end
        mean_rmses(strategy,cs_idx) = mean(rmses);
        std_rmses(strategy,cs_idx) = std(rmses);
        mean_lls(strategy,cs_idx) = mean(lls);
        std_lls(strategy,cs_idx) = std(lls);
        mean_ccs(strategy,cs_idx) = mean(ccs);
        std_ccs(strategy,cs_idx) = std(ccs);
    end
end
fprintf([fix_name(scen), ' & ']);
idx = 1;
n_output(mean_rmses(:,idx), 1);
fprintf(' & ');
n_output(mean_lls(:,idx), 0);
% for i=1:size(mean_rmses,1)
%     fprintf(strcat([' & %.2g', num2str(mean_rmses(i,end))]));
% end
fprintf('\\\\\n');
% mean_lls(:,end)'
 
return

figure
hold off
errorbar(caps,mean_rmses(3,:),std_rmses(3,:),std_rmses(3,:), 'gd', 'LineWidth',1, 'MarkerSize',10);
hold on
errorbar(caps,mean_rmses(4,:),std_rmses(4,:),std_rmses(4,:), 'r.', 'LineWidth',2, 'MarkerSize',20);
errorbar(caps,mean_rmses(2,:),std_rmses(2,:),std_rmses(2,:), 'bo', 'LineWidth',1, 'MarkerSize',10);
errorbar(caps,mean_rmses(1,:),std_rmses(1,:),std_rmses(1,:), 'kx', 'LineWidth',2, 'MarkerSize',10);
errorbarlogx
hx = xlabel(xLabelString);
hy = ylabel('RMSE');
legend({'drop censored', 'treat as uncensored', 'Schmee & Hahn', 'Sampling S&H'}, 'Location', 'NorthEast')
axis([caps(1)/1.3, caps(end)*16, min(min(mean_rmses-std_rmses))-0.1, max(max(mean_rmses+std_rmses))+0.1]);
set([gca], 'FontSize', 16);
set([hx,hy], 'FontSize', 20);
%set(gcf, 'Outerposition', [100,700,1000,480]);
set(gcf, 'Outerposition', [100,700,600,320]);
filename_prefix = strcat(['perf_capslack_rmse-', scen]);
saveas(gcf, strcat(filename_prefix,'.fig'));
set(gcf, 'PaperPositionMode', 'auto');
pdf_filename = strcat(filename_prefix,'.pdf');
if exist(pdf_filename, 'file');
    delete(pdf_filename);
end
export_fig(pdf_filename);

figure
errorbar(caps,mean_lls(1,:),std_lls(1,:),std_lls(1,:), 'kx', 'LineWidth',2, 'MarkerSize',10);
hold on
errorbar(caps,mean_lls(2,:),std_lls(2,:),std_lls(2,:), 'bo', 'LineWidth',1, 'MarkerSize',10);
errorbar(caps,mean_lls(4,:),std_lls(4,:),std_lls(4,:), 'r.', 'LineWidth',2, 'MarkerSize',20);
errorbar(caps,mean_lls(3,:),std_lls(3,:),std_lls(3,:), 'gd', 'LineWidth',1, 'MarkerSize',10);
errorbarlogx
hx=xlabel(xLabelString);
hy=ylabel('Log likelihood');
%legend({'Sampling S&H', 'Schmee & Hahn', 'treat as uncensored', 'drop censored'}, 'Location', 'SouthEast')
legend({'Sampling S&H', 'Schmee & Hahn', 'treat as uncensored', 'drop censored'}, 'Location', 'SouthEast')
axis([caps(1)/1.3, caps(end)*16, min(min(mean_lls-std_lls))-0.1, max(max(mean_lls+std_lls))+0.1]);
set([gca], 'FontSize', 16);
set([hx,hy], 'FontSize', 20);
%set(gcf, 'Outerposition', [0,0,1200,300]);
%set(gcf, 'Outerposition', [100,700,1000,400]);
set(gcf, 'Outerposition', [100,700,600,320]);

filename_prefix = strcat(['perf_capslack_ll-', scen]);
saveas(gcf, strcat(filename_prefix,'.fig'));
set(gcf, 'PaperPositionMode', 'auto');
pdf_filename = strcat(filename_prefix,'.pdf');
if exist(pdf_filename, 'file');
    delete(pdf_filename);
end
export_fig(pdf_filename);

% figure
% errorbar(caps,mean_ccs(1,:),std_ccs(1,:),std_ccs(1,:), 'kx');
% hold on
% errorbar(caps,mean_ccs(2,:),std_ccs(2,:),std_ccs(2,:), 'bd');
% errorbar(caps,mean_ccs(3,:),std_ccs(3,:),std_ccs(3,:), 'go');
% errorbar(caps,mean_ccs(4,:),std_ccs(4,:),std_ccs(4,:), 'r.');
% errorbarlogx
