function Mfun = build_preconditioner(A, precCfg)

    switch precCfg.type
        case 'ichol'
            L = ichol(A, precCfg.opts);
            Mfun = @(x) L \ (L' \ x);
        case 'none'
            Mfun = @(x) x;
        otherwise
        error('Unrecognized preconditioner type.');
    end

end
