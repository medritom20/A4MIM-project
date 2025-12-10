function dirPath = ensure_dir(dirPath)
    if ~exist(dirPath, 'dir')
        mkdir(dirPath);
    end
end