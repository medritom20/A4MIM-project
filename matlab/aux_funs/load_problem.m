function [A, problem] = load_problem(mtrxName)
    
    filePath = fullfile('./data',mtrxName) ;

    S = load(filePath);

    A = S.Problem.A;
    problem = S.Problem;

end
