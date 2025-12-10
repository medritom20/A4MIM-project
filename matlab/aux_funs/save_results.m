function save_results( matrixName, rngSeed, precCfg, results, sessions, testId)
    outputDir = ensure_dir( fullfile('.\out', 'results') );
    % Load existing results table if it exists
    resultsTableFile = fullfile(outputDir, sprintf('test%d_results_table.mat', testId));
    if isfile(resultsTableFile)
        load(resultsTableFile, 'resultsTable');
    else
        resultsTable = table();
    end

    %% UPDATE RESULTS TABLE
    % Append current results to table
    newRow = struct();
    newRow.timestamp = {datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss')};
    newRow.matrixName = {matrixName};
    newRow.rngSeed = rngSeed;
    newRow.precType = {precCfg.type};
    newRow.results = {results};
    newRow.sessions = {sessions};

    % Append computed results to results table
    resultsTable = [resultsTable; struct2table(newRow)];
    % Save updated results table
    save(resultsTableFile, 'resultsTable');
    fprintf('Results table updated: %s\n', resultsTableFile);

    %% SAVE INDIVIDUAL RESULTS FILE
    % Save individual results file
    timestamp = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
    precType = strrep(precCfg.type, '-', '_');
    fileName = sprintf('test%d_results_%s_seed%d_prec-%s_%s.mat', testId, matrixName(1:end-4), rngSeed, precType, timestamp);
    filePath = fullfile(outputDir, fileName);
    
    % Create structure with all results
    resultsData = struct();
    resultsData.matrixName = matrixName;
    resultsData.rngSeed = rngSeed;
    resultsData.precCfg = precCfg;
    resultsData.results = results;
    resultsData.sessions = sessions;
    
    save(filePath, 'resultsData');
    fprintf('Results data saved individually to: %s\n', filePath);
end