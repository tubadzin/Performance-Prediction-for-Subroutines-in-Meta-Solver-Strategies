function [res] = predict_submatrix(matrix, fullpredmatrix, fullpredmatrixvar, ind1, ind2, tuningScenario, suffix, modelType, learnTime, plotMatrices)
truematrix = matrix(ind2, ind1);
thispredmatrix = fullpredmatrix(ind2, ind1);
thispredmatrixvar = fullpredmatrixvar(ind2, ind1);

[rmse, ll, cc, cc_rank] = measures_of_fit(log10(truematrix(:)), log10(thispredmatrix(:)), thispredmatrixvar(:));
fprintf(strcat([modelType, ' performance for ', suffix, ': RMSE = %f, LL=%f, CC=%f, learnTime=%f\n']), [rmse, ll, cc, learnTime]); % CC=%f, cc, 
res = [rmse, ll, cc]; %, cc_rank];

if plotMatrices
    %=== For plotting, force predictions into the same
    %=== upper/lower limits.
    orig_thispredmatrix = thispredmatrix;
    thispredmatrix = max(thispredmatrix, 0.005);
    thispredmatrix = min(thispredmatrix, 300);
    
    figure
    h=image(log10(thispredmatrix), 'CDataMapping', 'scaled'); 
    colormap(gray); 
    hold on; 
    %title('mean pred matrix'); 
    if 0
        ylabel('config (sorted by true goodness)', 'FontSize', 22); 
        xlabel('instance (sorted by true hardness)', 'FontSize', 22);
    else
        if length(ind1) < 100
            set(gca, 'XTick', [1,length(ind1)], 'XTickLabel', {'easy','hard'})
        else
            set(gca, 'XTick', [50,length(ind1)-50], 'XTickLabel', {'easy','hard'})
        end
        if length(ind2) < 100
            set(gca, 'YTick', [1,length(ind2)], 'YTickLabel', {'good','bad'})
        else
            set(gca, 'YTick', [50,length(ind2)-50], 'YTickLabel', {'good','bad'})
        end
        
        ylabel('configurations', 'FontSize', 22); 
        xlabel('instances', 'FontSize', 22);
    end

    h=colorbar;
    mini = min(min(log10(thispredmatrix)));
    maxi = max(max(log10(thispredmatrix)));
    set(gca, 'FontSize', 16);
    set(h, 'YLim', [mini-1e-6,maxi+1e-6])
    mkdir('results/matrix/2d/');
    filename_prefix = strcat(['results/matrix/2d/', tuningScenario, '-predmatrix', suffix, '-', modelType]);
    saveas(gcf, strcat(filename_prefix,'.fig'));
    set(gcf, 'PaperPositionMode', 'auto');
    pdf_filename = strcat(filename_prefix,'.pdf');
    if exist(pdf_filename, 'file');
        delete(pdf_filename);
    end
    export_fig(pdf_filename);
    fprintf(strcat('Saved plot to ', pdf_filename, '\n'));

    figure
    h=image(log10(truematrix), 'CDataMapping', 'scaled'); 
    colormap(gray); 
    hold on; 
    %    title(tuningScenario); 
    if 0
        ylabel('config (sorted by true goodness)', 'FontSize', 22); 
        xlabel('instance (sorted by true hardness)', 'FontSize', 22);
    else
        if length(ind1) < 100
            set(gca, 'XTick', [1,length(ind1)], 'XTickLabel', {'easy','hard'})
        else
            set(gca, 'XTick', [50,length(ind1)-50], 'XTickLabel', {'easy','hard'})
        end
        if length(ind2) < 100
            set(gca, 'YTick', [1,length(ind2)], 'YTickLabel', {'good','bad'})
        else
            set(gca, 'YTick', [50,length(ind2)-50], 'YTickLabel', {'good','bad'})
        end
        ylabel('configurations', 'FontSize', 22); 
        xlabel('instances', 'FontSize', 22);
    end
    h=colorbar;
    mini = min(min(log10(truematrix)));
    maxi = max(max(log10(truematrix)));
    set(gca, 'FontSize', 16);
    set(h, 'YLim', [mini-1e-6,maxi+1e-6])
    filename_prefix = strcat(['results/matrix/2d/', tuningScenario, '-truematrix', suffix]);
    saveas(gcf, strcat(filename_prefix,'.fig'));
    set(gcf, 'PaperPositionMode', 'auto');
    pdf_filename = strcat(filename_prefix,'.pdf');
    if exist(pdf_filename, 'file');
        delete(pdf_filename);
    end
    export_fig(pdf_filename);

    %=== Plot part of the results in simple scatter plot.
    y = truematrix(:);
    ypred = orig_thispredmatrix(:);
    ypredvar = thispredmatrixvar(:);

    mkdir('results/matrix/single_preds/');
    filename_prefix = strcat(['results/matrix/single_preds/', tuningScenario, '-single_preds', suffix, '-', modelType]);
    perm = randperm(length(y));
    perm = perm(1:min(1000, length(perm)));
    title_prefix = tuningScenario;
    plot_simple_pred_scatter(y(perm), ypred(perm), ypredvar(perm), zeros(length(perm),1), rmse, cc, ll, filename_prefix, title_prefix, 0);

    % figure
    % trueMarginalMean = mean(truematrix, 2);
    % predMarginalMean = mean(thispredmatrix, 2);
    % plot(log10(trueMarginalMean), log10(predMarginalMean), '.');
    % 
    % [rmse, ll, cc, cc_rank] = measures_of_fit(log10(trueMarginalMean(:)), log10(predMarginalMean(:)), zeros(size(predMarginalMean(:),1),1));
    % rmse
end
