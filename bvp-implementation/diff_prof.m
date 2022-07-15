function [dfdx] = diff_prof(f, x)
    dfdx = 0.0*f;
    dfdx(2:(end-1)) = (f(3:end)-f(1:(end-2)))./(x(3:end)-x(1:(end-2)));
    dfdx(1) = (f(2)-f(1))/(x(2)-x(1));
    dfdx(end) = (f(end)-f(end-1))/(x(end)-x(end-1));
end

