function [A, problem] = load_problem(filePath)
    S = load(filePath);

    A = S.Problem.A;
    problem = S.Problem;

end
