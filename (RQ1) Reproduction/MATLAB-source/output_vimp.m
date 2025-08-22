function output_vimp(vimp, prefix)
for i=1:length(vimp.importance)
%         figure;
    [importance, idx] = sort(-vimp.importance{i});
    importance=-importance;
    feature = vimp.featureNames(vimp.selected_feature_idxs{i}(idx));
    
    figure;
    barh([importance]);
    set(gca,'YTickLabel',feature, 'FontSize', 14);
    title (sprintf('Relative Importance of Variables in Subset of Size %d (validation data)',size(importance,2)), 'FontSize', 14);
    xlabel ('Relative Importance', 'FontSize', 14);
    
    figure_prefix = strcat(prefix, '-', num2str(i));
    if ~strcmp(figure_prefix, '')
        set(gca, 'XTick', [0,20,40,60,80,100])
        axis([0 100 0.4 10.6])
%        set(gcf, 'y', 'auto');

        filename = strcat(figure_prefix, '-ehm-varImp.pdf');
        if exist(filename, 'file')
            delete(filename)
        end
        export_fig(filename);
        saveas(gcf, strcat(figure_prefix,'-ehm-varImp.fig'));
    end
    
%     figure
%     rmses = vimp.rmses{i};
%     h=plot(1:length(rmses), rmses);
%     hXLabel=xlabel ('Feature subset size', 'FontSize', 14);
%     hYLabel=ylabel ('RMSE', 'FontSize', 14);
% 
%     set(h, 'Color' , [1 1 1]);
%     set([hXLabel, hYLabel], 'FontSize', 16);
% 
%     set(h                      , ...
%     'LineStyle'       , ':'      , ...
%     'LineWidth'       , 2           , ...
%     'Color'           , [.4 .4 .4]  , ...
%     'MarkerSize'      , 6           , ...
%     'Marker'          , 'o'         , ...
%     'MarkerFaceColor' , [.7 .7 .7]  , ...
%     'MarkerEdgeColor' , [.4 .4 .4]);
% 
%     maxrmse = -inf;
%     minrmse = inf;
%     for i=1:length(vimp.rmses)
%         maxrmse = max(maxrmse, max(vimp.rmses{i}));
%         minrmse = min(minrmse, min(vimp.rmses{i}));
%         axis([0.001 11-0.001 0 maxrmse])
%     end

%         myboxplot( vimp.mean_importance(1:k), 'orientation', 'horizontal');
%         bar( vimp.mean_importance(1:k));
%         set(gca, 'XTick', 1:length(vimp.mean_importance(1:k)));
%         set(gca,'XTickLabel',vimp.feature(1:k), 'FontSize', 14);
%         
% 
% %        xlabel ('Relative Importance', 'FontSize', 14);
% 
%         figure_prefix = strcat('results/ehm/', domain, '-', num2str(i));
%         if ~strcmp(figure_prefix, '')
%             set(gcf, 'PaperPositionMode', 'auto');
% 
%             filename = strcat(figure_prefix, '-varImp.pdf');
%             if exist(filename, 'file')
%                 delete(filename)
%             end
%             export_fig(filename);
%             saveas(gcf, strcat(figure_prefix,'-varImp.fig'));
%         %    close;
%         end

end 

figure
hXLabel=xlabel ('Feature subset size', 'FontSize', 14);
hold on;
hYLabel=ylabel ('RMSE', 'FontSize', 14);

%set(h, 'Color' , [1 1 1]);
set([hXLabel, hYLabel], 'FontSize', 16);

for i=1:length(vimp.rmses)
    rmses = vimp.rmses{i};
    hE(i)=plot(1:length(rmses), rmses);
    hold on;
    hE2(i)=line([0,100],[vimp.full_rmse(i),vimp.full_rmse(i)]);
end

% set([hE(1)]                      , ...
set([hE(1),hE2(1)]                      , ...
  'LineStyle'       , ':'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [.4 .4 .4]  , ...
  'MarkerSize'      , 6           , ...
  'Marker'          , 'o'         , ...
  'MarkerFaceColor' , [.7 .7 .7]  , ...
  'MarkerEdgeColor' , [.4 .4 .4]);
set(hE2(1), 'LineStyle', '-');


if length(hE) > 1
%     set([hE(2)]                            , ...
    set([hE(2),hE2(2)]                            , ...
        'Color'           , 'b'  , ...
        'Marker'          , 'x'         , ...
        'LineStyle'       , '--'      , ...
        'MarkerSize'      , 12, ...        
        'LineWidth'       , 2);
    set(hE2(2), 'LineStyle', '-');

    if length(hE) > 2
%         set([hE(3)]                            , ...
        set([hE(3),hE2(3)]                            , ...
            'Color'           , 'r'  , ...
            'Marker'          , '.'         , ...
            'LineStyle'       , '-.'      , ...
            'MarkerSize'      , 16, ...        
            'LineWidth'       , 2);
%         set(hE2(3), 'LineStyle', '-');
        
        if length(hE) > 3
            set([hE(4),hE2(4)]                            , ...
                'Color'           , 'k'  , ...
                'Marker'          , 'd'         , ...
                'LineStyle'       , '-.'      , ...
                'MarkerSize'      , 9, ...        
                'LineWidth'       , 2);
            set(hE2(4), 'LineStyle', '-');

            if length(hE) > 4
                set([hE(5),hE2(5)]                            , ...
                    'Color'           , 'g'  , ...
                    'Marker'          , 's'         , ...
                    'LineStyle'       , ':'      , ...
                    'MarkerSize'      , 10, ...        
                    'LineWidth'       , 2);
                set(hE2(5), 'LineStyle', '-');
                
            end
        end
    end
end

maxrmse = -inf;
minrmse = inf;
legendstr = {};
for i=1:length(vimp.rmses)
    maxrmse = max(maxrmse, max(vimp.rmses{i}));
    minrmse = min(minrmse, min(vimp.rmses{i}));
    maxrmse = max(maxrmse, vimp.full_rmse(i));
    minrmse = min(minrmse, vimp.full_rmse(i));
    
    if isfield(vimp, 'options_vec')
        legendstr{i} = vimp.options_vec{i}.modelType;
    else
        legendstr{i} = 'a';
    end
end
%axis([0.001 11-0.001 max(0,minrmse-0.05) (maxrmse+0.05)*1.05])
axis([0.001 11-0.001 max(0,(minrmse-0.02)/1.01) (maxrmse*1.01)])

legendPos = 'NorthEast';
legend(hE, legendstr, 'location', legendPos);

figure_prefix = strcat(prefix, '-rmses_with_subsets');
if ~strcmp(figure_prefix, '')
%        set(gcf, 'y', 'auto');

    filename = strcat(figure_prefix, '-varImp.pdf');
    if exist(filename, 'file')
        delete(filename)
    end
    export_fig(filename);
    saveas(gcf, strcat(figure_prefix,'-varImp.fig'));
end

