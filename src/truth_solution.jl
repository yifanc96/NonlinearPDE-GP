alpha = 0;
m = 3;

# ground truth solution
function u(x1, x2)
    return sin(pi*x1)*sin(pi*x2) + 2*sin(4*pi*x1)*sin(4*pi*x2)
end

# right hand side
function f(x1, x2)
    return 2*pi^2*sin(pi*x1)*sin(pi*x2)+64*pi^2*sin(4*pi*x1)*sin(4*pi*x2)+alpha*u(x1,x2)^m
end

# boundary value
function g(x1, x2)
    return u(x1, x2)
end