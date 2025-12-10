function [mant, expo] = format_scientific(val)
% FORMAT_SCIENTIFIC  Returns mantissa and exponent for scientific notation.
    if val == 0
        mant = 0;
        expo = 0;
        return;
    end
    expo = floor(log10(abs(val)));
    mant = val / 10^expo;
    if abs(mant) >= 10
        mant = mant / 10;
        expo = expo + 1;
    elseif abs(mant) < 1
        mant = mant * 10;
        expo = expo - 1;
    end
end