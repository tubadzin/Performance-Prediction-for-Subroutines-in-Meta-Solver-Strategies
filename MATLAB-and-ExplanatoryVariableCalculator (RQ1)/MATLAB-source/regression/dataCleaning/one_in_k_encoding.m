function [encodedX, namesEncodedX] = one_in_k_encoding(X, cat, catDomains, namesX)
% Encode categorical parameters using a 1-in-k encoding.
if nargin < 4
    namesX = {};
end
cont = setdiff( 1:size(X,2), cat );
encodedX = X(:,cont);
namesEncodedX = {};
if ~isempty(namesX)
    namesEncodedX = namesX(cont);
end
for i=1:length(cat)
	x = X(:,cat(i));
	singleEncodedX = zeros(length(x), length(catDomains{cat(i)}));
	for j=1:length(x)
		singleEncodedX(j, x(j)) = 1;
    end
	encodedX = [encodedX, singleEncodedX];

    if ~isempty(namesX)
        singleEncodedNames = {};
        for j=1:length(catDomains{cat(i)})
            singleEncodedNames{end+1,1} = strcat(namesX{cat(i)}, '=', num2str(j));
        end
        namesEncodedX(end+1:end+length(catDomains{cat(i)}),1) = singleEncodedNames;
    end
end