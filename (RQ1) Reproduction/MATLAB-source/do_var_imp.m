function do_var_imp(options_vec, savefilename, X, y, cat, catDomains, featureNames)
k=10;
for i=1:k
    fprintf(strcat(['Sample for importance of input ', num2str(i), '/', num2str(k), '...\n']));
    this_savefilename = strcat(savefilename, '-', num2str(i));
    rand('twister', i);
    do_single_pass_of_var_imp(options_vec, this_savefilename, X, y, cat, catDomains, featureNames)
end