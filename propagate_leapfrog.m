function znext = propagate_leapfrog(f, x, xd, dt)
    xnext = x + xd*dt+0.5*f(x)*dt^2;
    xdnext = xd + 0.5*(f(x)+f(xnext))*dt;
    znext = [xnext;xdnext];
end