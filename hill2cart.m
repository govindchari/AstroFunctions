function [ri,vi] = hill2cart(rh,vh,t)
    a = 6778.137;
    mu = 3.986004418e5;
    n = sqrt(mu/a^3);
    w = [0;0;n];
    th = n*t;
    hCi = [cos(th) sin(th) 0;-sin(th) cos(th) 0;0 0 1];
    rOF_h = [a*cos(th);a*sin(th);0];
    rPF_h = rOF_h + rh;
    vPF_h = vh + cross(w,rPF_h);

    ri = hCi'*rPF_h;
    vi = hCi'*vPF_h;
end