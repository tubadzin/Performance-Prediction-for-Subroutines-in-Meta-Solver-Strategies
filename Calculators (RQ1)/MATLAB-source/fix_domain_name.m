function domain = fix_domain_name(domain)
domain = regexprep(domain, 'data_', '');
domain = regexprep(domain, '_csv', '');
domain = regexprep(domain, 'SAT_Competition_RACE_', '');
domain = regexprep(domain, 'SAT_', '');
domain = regexprep(domain, 'MIP_', '');
domain = regexprep(domain, 'TSP_', '');
domain = regexprep(domain, 'RAND_SAT', 'RAND-SAT');
