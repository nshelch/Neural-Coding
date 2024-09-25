function ai = computeAdaptationIndex(offFR, susFR, method)

if strcmpi(method, 'tan')
    ai = abs(atan(offFR ./ susFR) * (2/pi));
elseif strcmpi(method, 'norm')
    ai = offFR ./ (offFR + susFR);
else
    warning('Unsupported method requested')
end

end