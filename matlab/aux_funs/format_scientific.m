function txt = format_scientific(val)
% FORMAT_SCIENTIFIC  Returns a LaTeX-friendly string for scientific notation.

if val == 0
    txt = '0';
    return;
end

expo = floor(log10(abs(val)));
mant = val / 10^expo;

if abs(mant - 1) < 1e-12
    txt = sprintf('10^{%d}', expo);
else
    txt = sprintf('%g\\times 10^{%d}', mant, expo);
end

end
