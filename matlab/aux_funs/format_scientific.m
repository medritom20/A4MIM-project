function [mant, expo] = format_scientific(value)
% FORMAT_SCIENTIFIC  Returns mantissa and exponent for scientific notation.

    if value < 0
        error('format_scientific:NegativeValue', 'Input value must be non-negative.');
    elseif value == 0
        mant = 0;
        expo = 0;
        return;
    end

    expo = floor(log10(abs(value)));
    mant = value / 10^expo;
end