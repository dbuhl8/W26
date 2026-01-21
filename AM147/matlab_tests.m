% Newtons Method Test 2

h = 1e-5;
f= @(x) 2*x^x - 7*log(x) - 3;
fprime = @(x,f) (f(x+h) - f(x-h))/(2*h);
y = NewtonsMethod(f,fprime,1,1e-8);


% Newtons Method Test 3 (looks like a typo hehe)
h = 1e-5;
f= @(x) e^x^2 - .001*log(x) - 10;% Run learner solution.
h = 1e-5;
f= @(x) 2*x^x - 7*log(x) - 3;
fprime= @(x,f) (f(x+h) - f(x-h))/(2*h);
y = NewtonsMethod(f,fprime,1,1e-8);
