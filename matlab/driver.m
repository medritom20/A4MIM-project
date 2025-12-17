function results = driver(sessions)
% DRIVER  Orchestrates solver runs for a list of sessions.
%
% Each element of SESSIONS describes a block-size experiment:
%   .matrixName    matrix name under ./data
%   .blockSize     block size m (#RHSs)
%   .rngSeed       (optional) RNG seed for reproducibility
%   .prec          preconditioner configuration (type/opts)
%   .runs          array with fields:
%                     method
%                     policy  - solver parameters (tol, maxIt, etc.)
%                     label   - optional legend label override
%                     style   - optional plot style override
%
% The function returns a RESULTS struct containing the per-session outputs.

    [A, ~] = load_problem(sessions(1).matrixName);
    n = size(A, 1);
    [applyM, precInfo] = build_preconditioner(A, sessions(1).prec);

    % results.problem = problem;
    results.precInfo = precInfo;
    nSessions = numel(sessions);
    % results.sessions = repmat(struct(), numSessions, 1);

    for iS = 1:nSessions
        sess = sessions(iS);

        m = sess.blockSize;
        rngSeed = sess.rngSeed + iS - 1;
        rng(rngSeed, 'twister');

        B     = randn(n, m); 
        Xtrue = A \ B;      % reference solution via direct solve 
        X0    = zeros(n, m);

        nRuns = numel(sess.runs);
        fprintf('\nSession %2d/%2d (m = %d, rngSeed = %d)\n', iS, nSessions, m, rngSeed);

        % runResults = runsCfg;
        for im = 1:nRuns
            runCfg = sess.runs(im);
            methodName = sess.runs(im).method;
            policy = sess.runs(im).policy;
            fprintf('\tRun %d/%d: %s\n', im, nRuns, methodName);
            runOutput = run_method(methodName, A, B, X0, Xtrue, applyM, precInfo, policy);

            runEntry        = runCfg;
            runEntry.name   = methodName;
            runEntry.method = methodName;
            runEntry.policy = policy;
            runEntry.results = runOutput;
            runResults(im)  = runEntry;
        end

        sess.runs = reshape(runResults, size(sess.runs));
        sess.blockSize = m;
        sess.rngSeed = rngSeed;
        results.sessions(iS) = sess;
    end

end




function runOut = run_method(name, A, B, X0, Xtrue, applyM, precInfo, policy)
    switch name
        case 'HS-BCG'
            if ~isfield(policy,'phiType') || isempty(policy.phiType)
                policy.phiType = 'eye';
            end
            [~, flag, relres, iter, resHist, omegaHist] = HS_PBCG( ...
                A, B, applyM, X0, policy.tol, policy.maxIt, policy.phiType, Xtrue);
            runOut = struct('omega', omegaHist, 'flag', flag, 'relres', relres, ...
                            'iter', iter, 'resHist', resHist);

        case 'DP-BCG'
            [~, omegaHist] = PDPBCG(A, B, X0, applyM, Xtrue, policy.maxIt, policy.tol);
            runOut = struct('omega', omegaHist);

        case 'DR-BCG'
            if strcmpi(precInfo.type, 'none')
                [omegaHist, ~] = DR_BCG_exp(A, B, X0, Xtrue, policy.maxIt);
            else
                if ~isfield(precInfo, 'L') || isempty(precInfo.L)
                    error('driver:missingL', 'Preconditioned DR-BCG requires the Cholesky factor L.');
                end
                [omegaHist, ~] = Prec_DR_BCG_exp(A, B, X0, Xtrue, precInfo.L, policy.maxIt);
            end
            runOut = struct('omega', omegaHist);

        case 'BF-BCG'
            if ~isfield(policy, 'svdTol')
                error('driver:missingSvdTol', 'Missing svdTol parameter for BF-BCG.');
            end
            % [~, omegaHist] = BFBCG(A, B, X0, applyM, Xtrue, policy.maxIt, policy.tol, svdTol);
            % [~, ~, ~, omegaHist] = pbf_bcg(A, B, applyM, Xtrue, policy.maxIt, policy.tol, svdTol);
            [~, omegaHist] = PBFBCG(A, B, X0, applyM, Xtrue, policy.maxIt, policy.tol, policy.svdTol);
            runOut = struct('omega', omegaHist);

        otherwise
            error('driver:unknownMethod', 'Unknown method "%s".', name);
    end
end
