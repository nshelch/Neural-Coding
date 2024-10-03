function y = gamma_fx(abc, x)
a = abc(1); b = abc(2); c = abc(3);
y = c * x.^(a-1) .* exp(-x/b) / (b^a * gamma(a));
end