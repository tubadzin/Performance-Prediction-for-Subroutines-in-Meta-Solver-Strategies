function output_vimps(savefile, prefix, maxPlotModels, xstring)
if nargin < 4
    xstring = 'Feature subset size';
    if nargin < 3
        maxPlotModels = inf;
    end
end
vimps = {};
numSamples = 10;
for n=1:numSamples
    load(strcat([savefile, '-' num2str(n)]), 'vimp');
    vimps{n} = vimp;
end
len = length(vimps{1}.options_vec);
len = min(len, maxPlotModels);

featureNames = vimps{1}.featureNames;
numFeatures = length(featureNames);
for model_idx = 1:len
    importance = zeros(numFeatures, numSamples);
    for n = 1:numSamples
        importance(vimps{n}.selected_feature_idxs{model_idx}, n) = vimps{n}.importance{model_idx};
        importance(find(importance(:,n)>99.99),n) = 100;
        importance(find(importance(:,n)<0.01),n) = 0;
    end
    mean_importance = mean(importance,2);
    [tmp, sort_idx] = sort(-mean_importance);
    k = 10;
    idx = sort_idx(k:-1:1);
    figure
    myboxplot( importance(idx(1:k),:)', 'orientation', 'horizontal' );
%     for i=1:k
%         plot(importance(idx(i),:), ones(numSamples,1)*i, 'ko');
%         hold on
%     end
	set(gca, 'YTick', 1:k);
	set(gca,'YTickLabel',featureNames(idx(1:k)), 'FontSize', 14);
    ylabel('');
    xlabel('Importance');
    axis([-5 105 0.001 k+1-0.001])
    
    figure_prefix = strcat(prefix, '-', num2str(model_idx));
    filename = strcat(figure_prefix, '-ehm-varImp.pdf');
    if exist(filename, 'file')
        delete(filename)
    end
    export_fig(filename);
    saveas(gcf, strcat(figure_prefix,'-ehm-varImp.fig'));
end


figure
hXLabel=xlabel(xstring);
hold on;
hYLabel=ylabel('RMSE');

plot_full_rmse = 1;

maxrmse = -inf;
minrmse = inf;
hE = [];
hE2 = [];
legendstr = {};
for model_idx = 1:len
    mean_rmses = zeros(1, length(vimps{1}.rmses{1}));
    full_rmse = 0;
    for numSample = 1:numSamples
        rmses = vimps{numSample}.rmses{model_idx};
        mean_rmses = mean_rmses + rmses;
        full_rmse = full_rmse + vimps{numSample}.full_rmse(model_idx);
    end
    mean_rmses = mean_rmses/numSamples;
    full_rmse = full_rmse/numSamples;
    hE(end+1)=plot(1:length(mean_rmses), mean_rmses);
    hold on;
    maxrmse = max(maxrmse, max(mean_rmses));
    minrmse = min(minrmse, min(mean_rmses));
    if plot_full_rmse
        hE2(end+1)=line([0,100],[full_rmse,full_rmse]);
        maxrmse = max(maxrmse, full_rmse);
        minrmse = min(minrmse, full_rmse);
    end
    if isfield(vimp, 'options_vec')
        legendstr{end+1} = vimp.options_vec{model_idx}.modelType;
        if strcmp(legendstr{end}, 'rf')
            legendstr{end} = 'Random forest';
        end
        if strcmp(legendstr{end}, 'GPML')
            legendstr{end} = 'Projected process';
        end
        if strcmp(legendstr{end}, 'LR')
            legendstr{end} = 'Ridge regression';
        end
    else
        legendstr{end+1} = 'a';
    end
end

% for i=1:length(vimp.rmses)
%     rmses = vimp.rmses{i};
%     hE(i)=plot(1:length(rmses), rmses);
%     hold on;
%     hE2(i)=line([0,100],[vimp.full_rmse(i),vimp.full_rmse(i)]);
% end

set(hE(1)                      , ...
  'LineStyle'       , '--'      , ...
  'LineWidth'       , 2           , ...
  'Color'           , [.4 .4 .4]  , ...
  'MarkerSize'      , 6           , ...
  'Marker'          , 'o'         , ...
  'MarkerFaceColor' , [.7 .7 .7]  , ...
  'MarkerEdgeColor' , [.4 .4 .4]);
if plot_full_rmse
    set(hE2(1), 'LineStyle', '-', ...
  'LineStyle'       , '--'      , ...
  'LineWidth'       , 1           , ...
  'Color'           , [.4 .4 .4]  , ...
  'MarkerSize'      , 6           , ...
  'Marker'          , 'o'         , ...
  'MarkerFaceColor' , [.7 .7 .7]  , ...
  'MarkerEdgeColor' , [.4 .4 .4]);
end


if length(hE) > 1
    set(hE(2)                            , ...
        'Color'           , 'b'  , ...
        'Marker'          , 'x'         , ...
        'LineStyle'       , '-'      , ...
        'MarkerSize'      , 12, ...        
        'LineWidth'       , 2);
    if plot_full_rmse
        set(hE2(2), 'LineStyle', '-', ...
        'Color'           , 'b'  , ...
            'Marker'          , 'x'         , ...
            'MarkerSize'      , 12, ...        
            'LineWidth'       , 1);
    end

    if length(hE) > 2
        set(hE(3)                            , ...
            'Color'           , 'r'  , ...
            'Marker'          , '.'         , ...
            'LineStyle'       , ':'      , ...
            'MarkerSize'      , 16, ...        
            'LineWidth'       , 2);
        if plot_full_rmse
             set(hE2(3), 'LineStyle', ':',  ...
            'Color'           , 'r'  , ...
            'Marker'          , '.'         , ...
            'MarkerSize'      , 16, ...        
            'LineWidth'       , 1);
        end
        
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

legendstr = legendstr(len:-1:1);
hE = hE(len:-1:1);

slack = (maxrmse-minrmse) / 10;
%axis([0.001 11-0.001 max(0,minrmse-0.05) (maxrmse+0.05)*1.05])
%axis([0.001 11-0.001 max(0,(minrmse-slack)/1) (maxrmse+slack)*1])
axis([0.001 11-0.001 0 (maxrmse+slack)*1])

legendPos = 'SouthWest';
legend(hE, legendstr, 'location', legendPos);

fontSize = 24;
set([hXLabel, hYLabel], 'FontSize', fontSize);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'FontSize'    , fontSize-10, ...
  'LineWidth'   , 2         );
set(gcf, 'Outerposition', [0,0,500,500]);
%   'XTick'       , 10.^[-10:10], ...
%   'YTick'       , 10.^[-10:10], ...

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

