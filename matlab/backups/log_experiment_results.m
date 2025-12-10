function resultsTable = log_experiment_results(logFilePath, rngSeed, policyCfg, results)
% LOG_EXPERIMENT_RESULTS  Persist a single experiment run into a MAT table.
%
% resultsTable = LOG_EXPERIMENT_RESULTS(logFilePath, rngSeed, policyCfg, results)
% creates/updates a MAT-file located at LOGFILEPATH. The file stores a table
% named RESULTS TABLE with the following columns:
%   rngSeed   Base RNG seed that produced the experiment.
%   policy    Copy of the session configuration (policies) used for the run.
%   results   Driver output, including per-session results.
%   timestamp Timestamp indicating when the entry was created.
%
% When an entry with the same rngSeed already exists it is overwritten with
% the newest data so repeated runs simply refresh their record.

    arguments
        logFilePath {mustBeTextScalar}
        rngSeed (1,1) double
        policyCfg
        results
    end

    logDir = fileparts(logFilePath);
    if isempty(logDir)
        logDir = '.';
    end
    if ~exist(logDir, 'dir')
        mkdir(logDir);
    end

    resultsTable = table();
    if isfile(logFilePath)
        loaded = load(logFilePath, 'resultsTable');
        if isfield(loaded, 'resultsTable')
            resultsTable = loaded.resultsTable;
        end
    end

    newRow = table( ...
        rngSeed, {policyCfg}, {results}, datetime('now'), ...
        'VariableNames', {'rngSeed', 'policy', 'results', 'timestamp'});

    if isempty(resultsTable)
        resultsTable = newRow;
    else
        idx = resultsTable.rngSeed == rngSeed;
        if any(idx)
            resultsTable(idx, :) = newRow;
        else
            resultsTable = [resultsTable; newRow]; %#ok<AGROW>
        end
    end

    save(logFilePath, 'resultsTable');
end

function mustBeTextScalar(value)
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        error('log_experiment_results:invalidPath', ...
              'logFilePath must be a character vector or string scalar.');
    end
end
