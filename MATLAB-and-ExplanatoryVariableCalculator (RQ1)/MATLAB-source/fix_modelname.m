function name = fix_modelname(name)
switch name
    case 'GP'
        name = 'PP';
    case 'RR-cv'
        name = 'RR(cv)';
    case 'RF-def'
        name = 'RF';
    case 'RT'
        name = 'RT';
    case 'RR'
        name = 'RR';
    case 'RR-drop-linear-comb'
        name = 'RR-elim';
    case 'Foba-cv-rmse'
        name = 'SP(cv)';
    case 'Foba-rmse'
        name = 'SP';
    case 'NN'
        name = 'NN';
    case 'NN-cv'
        name = 'NN(cv)';
    case 'RF-cv'
        name = 'RF(cv)';
end        