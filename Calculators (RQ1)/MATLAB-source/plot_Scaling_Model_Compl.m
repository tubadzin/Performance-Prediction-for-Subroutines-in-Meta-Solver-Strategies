function plot_Scaling_Model_Compl(numScal, namesNumScal, data_avg, data_std, modelTimesAvg, dataname, model_names, domain, legendPos, y_log)
if nargin < 10
    y_log = 0;
    if nargin < 9
        legendPos = 'SouthEast';
    end
end

% hE = [];
% for model_idx=1:length(data_avg)
%     tmp_avg = [];
%     tmp_std = [];
%     tmp_numscal = [];
%     for i=1:length(data_avg{model_idx})
%         tmp_avg(i) = data_avg{model_idx}{i};
%         tmp_std(i) = data_std{model_idx}{i};
%         tmp_numscal(i) = numScal{model_idx}(i);
%     end
%     
%     figure
%     hold off
%     if y_log
%         h = loglog(tmp_numscal, tmp_avg);
%     else
%         h = semilogx(tmp_numscal, tmp_avg);
%     end
%     
%     hold on
%     hE(end+1) = errorbar(tmp_numscal, tmp_avg, tmp_std, 'r');
%     hXLabel = xlabel(namesNumScal{model_idx});
%     hYLabel = ylabel(dataname);
% 
%     mintime = min(tmp_numscal)/1.3;
%     maxtime = max(tmp_numscal)*1.3;
%     
%     set(gca, ...
%     'XTick'       , tmp_numscal, ...
%     'FontName'   , 'Helvetica' , ...
%     'FontSize'   , 16);
% %    'YTick'       , -1:0.1:1   , ...
% 
%     max_data = max(max(tmp_avg+tmp_std + 0.03));
%     min_data = min(min(tmp_avg-tmp_std - 0.03));
% 
%     axis([mintime, maxtime, min_data, max_data])
%     
%     set(h, 'Color' , [1 1 1]);
%     set([hXLabel, hYLabel], 'FontSize', 16);
% 
%     set(hE(1)                      , ...
%       'LineStyle'       , ':'      , ...
%       'LineWidth'       , 2           , ...
%       'Color'           , [.4 .4 .4]  , ...
%       'MarkerSize'      , 6           , ...
%       'Marker'          , 'o'         , ...
%       'MarkerFaceColor' , [.7 .7 .7]  , ...
%       'MarkerEdgeColor' , [.4 .4 .4]);
% 
%     if length(hE) > 1
%         set([hE(2)]                            , ...
%             'Color'           , 'b'  , ...
%             'Marker'          , 'x'         , ...
%             'MarkerSize'      , 8, ...        
%             'LineWidth'       , 2);
% 
%         if length(hE) > 2
%             set([hE(3)]                            , ...
%                 'Color'           , 'r'  , ...
%                 'Marker'          , '.'         , ...
%                 'LineStyle'       , '--'      , ...
%                 'MarkerSize'      , 16, ...        
%                 'LineWidth'       , 2);
%         end
%     end
% 
%     legend(hE(model_idx), model_names(model_idx), 'location', legendPos);
%     filename = strcat('scale_Compl_', domain, '_', dataname, '_', model_names{model_idx});
%     pdf_filename = strcat(filename, '.pdf');
%     if exist(pdf_filename, 'file');
%         delete(pdf_filename);
%     end
%     export_fig(pdf_filename);    
%     saveas(gca, strcat(filename, '.fig'), 'fig');    
% end


%=== Rest of file: x axis is model time, all models combined in one plot

hE = [];
legend_entries = {};
max_data = -inf;
maxtime = -inf;
min_data = inf;
mintime = inf;
for model_idx=1:length(data_avg)
    tmp_avg = [];
    tmp_std = [];
    tmp_times = [];
    for i=1:length(data_avg{model_idx})
        tmp_avg(i) = data_avg{model_idx}(i);
        tmp_std(i) = data_std{model_idx}(i);
        tmp_times(i) = modelTimesAvg{model_idx}(i);
    end
    if model_idx == 1
        figure
        hold off
        if y_log
            h = loglog(tmp_times, tmp_avg);
        else
            h = semilogx(tmp_times, tmp_avg);
        end
        hold on
    end
    
    hE(end+1) = errorbar(tmp_times, tmp_avg, tmp_std, 'r');
    legend_entries{end+1} = model_names(model_idx);
    hXLabel = xlabel('Model time (learn+predict [s])');
    hYLabel = ylabel(dataname);
    
    max_data = max(max_data, max(max(tmp_avg+tmp_std + 0.03)));
    min_data = min(min_data, min(min(tmp_avg-tmp_std - 0.03)));
    
    mintime = min(mintime, min(tmp_times))/1.3;
    maxtime = max(maxtime, max(tmp_times))*1.3;
end

set(h, 'Color' , [1 1 1]);
set([hXLabel, hYLabel], 'FontSize', 16);

set(hE(1)                      , ...
  'LineStyle'       , '-'      , ...
  'LineWidth'       , 1           , ...
  'Color'           , 'k'  , ...
  'MarkerSize'      , 6           , ...
  'Marker'          , 'o'         , ...
  'MarkerFaceColor' , [.7 .7 .7]  , ...
  'MarkerEdgeColor' , 'k'); %[.4 .4 .4]

if length(hE) > 1
    set([hE(2)]                            , ...
        'Color'           , 'b'  , ...
        'Marker'          , 'x'         , ...
        'LineStyle'       , ':'      , ...
        'MarkerSize'      , 15, ...        
        'LineWidth'       , 2);

    if length(hE) > 2
        set([hE(3)]                            , ...
            'Color'           , 'r'  , ...
            'Marker'          , '+'         , ...
            'MarkerSize'      , 10, ...        
            'LineWidth'       , 2, ...
            'LineStyle', ':');
        
        if length(hE) > 3
            set([hE(4)]                            , ...
                'Color'           , 'g'  , ...
                'Marker'          , 'd'         , ...
                'LineStyle'       , '-'      , ...
                'MarkerSize'      , 9, ...        
                'LineWidth'       , 1);

%                    'Color'           , [.4 .4 .4]  , ...            
            if length(hE) > 4
                set([hE(5)]                            , ...
                    'Color'           , 'm'  , ...
                    'Marker'          , 's'         , ...
                    'LineStyle'       , '-'      , ...
                    'MarkerSize'      , 10, ...        
                    'LineWidth'       , 1);
            
                if length(hE) > 5
                    set([hE(6)]                            , ...
                        'Color'           , 'k'  , ...
                        'Marker'          , 'o'         , ...
                        'LineStyle'       , ':'      , ...
                        'MarkerSize'      , 10, ...        
                        'LineWidth'       , 1);

                    if length(hE) > 6
                        set([hE(7)]                            , ...
                            'Color'           , [.4 .4 .4]  , ...
                            'Marker'          , 'h'         , ...
                            'LineStyle'       , '-'      , ...
                            'MarkerSize'      , 10, ...        
                            'LineWidth'       , 2);
                    end
                end
            end            
        end
    end
end

set(gca, ...
    'XTick'       , 10.^[-5:5], ...
    'FontName'   , 'Helvetica' , ...
    'FontSize'   , 16);
%    'YTick'       , -1:0.1:1   , ...

axis([mintime, maxtime, min_data, max_data])

legend_str = legend_entries{1};
for i=2:length(legend_entries)
    legend_str = strcat([legend_str, legend_entries{i}]);
end
legend(hE, legend_str, 'location', legendPos);
dataname = regexprep(dataname, ' ', '_');
filename = strcat(domain, '_', dataname);
pdf_filename = strcat(filename, '.pdf');
if exist(pdf_filename, 'file');
    delete(pdf_filename);
end
export_fig(pdf_filename);    
saveas(gca, strcat(filename, '.fig'), 'fig');