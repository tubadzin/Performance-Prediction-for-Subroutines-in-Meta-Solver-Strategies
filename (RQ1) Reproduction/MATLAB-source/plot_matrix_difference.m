names = {};
names{end+1} = 'CPLEX12-cat-CORLAT';
names{end+1} = 'CPLEX12-cat-BIGMIX';
names{end+1} = 'CPLEX12-cat-REG';
names{end+1} = 'CPLEX12-cat-CORLAT-REG';
names{end+1} = 'CPLEX12-cat-RCW';
names{end+1} = 'CPLEX12-cat-CORLAT-REG-RCW';
names{end+1} = 'SPEAR-ibm-al';
names{end+1} = 'SPEAR-swv-al';
names{end+1} = 'SPEAR-ibm-swv-al';

elements = {'trainC_trainI', 'trainC_testI', 'testC_trainI', 'testC_testI'};

for j=1:length(names)
    name = names{j}
    for i=1:length(elements)
        el = elements{i};

    	fighandle=openfig(['/ubc/cs/project/arrow/hutter/EHM-Code/results/matrix/2d/',name,'-100000-predmatrix-', el, '-rf.fig']);
%         fighandle=openfig(['/ubc/cs/project/arrow/hutter/EHM-Code/results/matrix/2d/',name,'-predmatrix-', el, '-rf.fig']);
        fighandle2=openfig(['/ubc/cs/project/arrow/hutter/EHM-Code/results/matrix/2d/',name,'-truematrix-', el, '.fig']);

        a_pred=getimage(fighandle);
        a_true=getimage(fighandle2);
        % close all;

        % mmax = max(max(abs(a_true-a_pred)));
        % h=image(mmax-abs(a_true-a_pred), 'CDataMapping', 'scaled');
        % colormap(gray);
        figure
%         h=image(abs(a_true-a_pred), 'CDataMapping', 'scaled');
        matrix = a_true-a_pred;
        matrix(1,1) = -4.5;
        matrix(end,end) = +4.5;
        
        h=image(matrix, 'CDataMapping', 'scaled');
        colormap(hot);
        h=colorbar;

        set(gca, 'XTick', [50,size(a_true,2)-50], 'XTickLabel', {'easy','hard'})
        set(gca, 'YTick', [50,size(a_true,1)-50], 'YTickLabel', {'good','bad'})
        ylabel('configurations', 'FontSize', 22); 
        xlabel('instances', 'FontSize', 22);

        set(gca, 'FontSize', 16);
%         set(h, 'YLim', [0,4.5]);
%         set(h, 'YLim', [min(a_true-a_pred),max(a_true-a_pred)]);
%         set(h, 'YLim', [-4.5,4.5]);

        filename_prefix = strcat(['results/matrix/2d/redone/', name, '-', el]);
        saveas(gcf, strcat(filename_prefix,'.fig'));
        set(gcf, 'PaperPositionMode', 'auto');

        pdf_filename = strcat(filename_prefix,'.pdf');
        if exist(pdf_filename, 'file');
            delete(pdf_filename);
        end
        export_fig(pdf_filename);
        fprintf(strcat('Saved plot to ', pdf_filename, '\n'));
    end
end