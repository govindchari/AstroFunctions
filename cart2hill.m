function [rh,vh] = cart2hill(ri,vi,t)
    a = 6778.137;
    mu = 3.986004418e5;
    n = sqrt(mu/a^3);
    w = [0;0;n];
    th = n*t;
    hCi = [cos(th) sin(th) 0;-sin(th) cos(th) 0;0 0 1];
    rOF_h = [a*cos(th);a*sin(th);0];

    rh = hCi*ri - rOF_h;
    vh = hCi*vi - cross(w,hCi*ri);
end