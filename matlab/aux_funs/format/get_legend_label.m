function label = get_legend_label(prec)
    if ~isfield(prec, 'type') || strcmpi(prec.type, 'none')
        label = 'DP-BCG';
        return;
    end

    switch lower(prec.type)
        case 'ichol'
            dropTol = prec.opts.droptol;
            [mant, expo] = format_scientific(dropTol);
            if mant == 1
                tolStr = sprintf('10^{%d}', expo);
            else
                tolStr = sprintf('%g\\times10^{%d}', mant, expo);
            end
            label = sprintf('PDP-BCG : ICT($%s$)', tolStr);
        otherwise
            label = sprintf('\\text{PDP-BCG %s}', upper(prec.type));
    end
end