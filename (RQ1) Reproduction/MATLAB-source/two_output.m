function two_output(a,b,small_good)
a=two_dec(a);
b=two_dec(b);

if ((a < b) && (small_good==1)) || ((a > b) && (small_good==0))
    fprintf(strcat(['$\\bm{', to_s(a), '}$ & ', to_s(b)]));
elseif ((a > b) && (small_good==1)) || ((a < b) && (small_good==0))
    fprintf(strcat([to_s(a), ' & $\\bm{', to_s(b), '}$']));
else
    fprintf(strcat(['$\\bm{', to_s(a), '}$ & $\\bm{', to_s(b), '}$']));
end

function x=two_dec(x)
x = ceil(x*100-0.5)/100;

function x = to_s(x)
x = sprintf( '%4.2f' , x );