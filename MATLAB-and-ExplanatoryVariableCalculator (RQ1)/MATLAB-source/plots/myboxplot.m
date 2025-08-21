function h=myboxplot(varargin)
h = boxplot(varargin{:}, 'colors', 'kkkk', 'symbol', 'ko');

for i=1:size(h,1)
    for j=1:size(h,2)
        if ~isnan(h(i,j))
            set(h(i,j),'LineWidth',1)
        end
    end
end
