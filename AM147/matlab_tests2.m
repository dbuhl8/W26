% Bisection Method Test 2

f= @(x) log(3*x^2) + x/50;
y = BisectionMethod(f,0,2,1e-10)
% Run reference solution.
yReference = reference.BisectionMethod(f,0,2,1e-10); 


% Bisection Method Test 3
f= @(x) 3*x-9;
y = BisectionMethod(f,2,4,1e-10)
% Run reference solution.
yReference = reference.BisectionMethod(f,2,4,1e-10); 
