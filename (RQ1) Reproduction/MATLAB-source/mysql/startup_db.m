function startup_db(func)
if ~(func.cheap || func.matlab_fun)
    if mysql('status')
        startUpDB
    end
    if func.mip
        mysql('use MIP');
    else
        mysql('use UBCSAT')
%        mysql('use tmptest')
%        mysql('use newtmptest')
    end
end