function plot_simple_pred_scatter(y, y_cross, y_cross_var, cens, rmse, cc, ll, figure_prefix, title_prefix, logModel, axmin, axmax, fontSize)
if nargin < 13
    fontSize = 30;
end
nCross = length(y_cross);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% observed values vs. cross-validated prediction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cens_idx=find(cens==1);
uncens_idx=find(cens==0);

if logModel
    ycross_plot = 10.^y_cross;
else
    ycross_plot = y_cross;
end
if nargin < 12
            axmin = 0.001;
            axmax = 3600;
%     axmin = min(min(y), min(ycross_plot));
%     axmax = max(max(y), max(ycross_plot));
end


margin_on_each_side = 1.3;

sfigure;
yplot = y;
if logModel
    ycross_plot = 10.^y_cross;
    if max(ycross_plot) > axmax+1e-5
        fprintf(strcat(['\n\nWARNING\n\nAt least one entry larger than upper axis boundary:', num2str(max(ycross_plot)), '\n\n']));
        idx = find(ycross_plot > axmax);
        num_preds_above_upper_boundary = length(idx)
        fprintf(strcat(['\nEND OF WARNING.\n\n']));
    end
    if min(ycross_plot) < axmin-1e-5
        fprintf(strcat(['\n\nWARNING\n\nAt least one entry smaller than lower axis boundary:', num2str(min(ycross_plot)), '\n\n']));
        idx = find(ycross_plot < axmin-1e-5);
        num_preds_below_lower_boundary = length(idx)
        fprintf(strcat(['\nEND OF WARNING.\n\n']));
    end
    bad_idx = find(ycross_plot(:) < axmin-1e-5);
    bad_idx = [bad_idx; find(ycross_plot(:) > axmax)];
else
    ycross_plot = y_cross;
    bad_idx = find(ycross_plot(:) < axmin);
    bad_idx = [bad_idx; find(ycross_plot(:) > axmax)];
    ycross_plot = y_cross;
end
good_idx = setdiff(1:length(ycross_plot), bad_idx);
mini = axmin/margin_on_each_side;
maxi = axmax*margin_on_each_side;

hold off
hE     = loglog(yplot(good_idx), ycross_plot(good_idx));
hold on

hLine = line([mini, maxi],[mini,maxi]);

set(hLine                         , ...
  'Color'           , [0 0 .5]    , ...
  'LineWidth'       , 2           );

set([hE]                     , ...
  'LineStyle'       , 'none'      , ...
  'Color'           , [.3 .3 .3]  , ...
  'LineWidth'       , 1           , ...
  'MarkerSize'      , 3           , ...
  'Marker'          , 'o'         , ...
  'MarkerFaceColor' , [.7 .7 .7]  , ...
  'MarkerEdgeColor' , [.2 .2 .2]);

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XColor'      , [0 0 0], ...
  'YColor'      , [0 0 0], ...
  'FontSize'    , fontSize-4, ...
  'XTick'       , 10.^[-10:6], ...
  'YTick'       , 10.^[-10:6], ...
  'LineWidth'   , 2         );


if ~isempty(bad_idx)
    ycross_plot = min(ycross_plot, axmax);
    ycross_plot = max(ycross_plot, axmin);
    hE2     = loglog(yplot(bad_idx), ycross_plot(bad_idx), 'bx', 'MarkerSize', 6);
end

axis([mini, maxi, mini, maxi]);

set(gcf, 'Outerposition', [0,0,500,500]);

if ~strcmp(figure_prefix,'')
    filenamePDF = strcat(figure_prefix,'-pred.pdf');
    filenameFIG = strcat(figure_prefix,'-pred.fig');
    exportgraphics(gcf, filenamePDF, 'ContentType','vector');
    savefig(gcf, filenameFIG);
end