function p = gammaFit(x, y)

p = nlinfit(x, y, @gamma_fx, [1 10 150]); %[1 10 150]

end


