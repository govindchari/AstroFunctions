function [r,v] = kepler2rv(a,e,Omega,I,omega,tmtp,mu)
    %INPUTS
    %   a - semjax 
    %   e - eccentricity
    %   Omega - longitude of ascending node
    %   I - inclination
    %   omega - argument of periapsis
    %   tmtp - time since periapsis passage (t-tp)
    %   mu - gravitational parameter (units of this dictate units of output)
    %OUTPUTS
    %   r - inertial position vector
    %   v - inertial velocity vector

	b = a*sqrt(1-e^2);
	CO = [cos(Omega) sin(Omega) 0;-sin(Omega) cos(Omega) 0;0 0 1];
	CI = [1 0 0;0 cos(I) sin(I);0 -sin(I) cos(I)];
	Co = [cos(omega) sin(omega) 0;-sin(omega) cos(omega) 0;0 0 1];
	pCi = Co*CI*CO;
	n = sqrt(mu/a^3);
	M = mod(n*(tmtp),2*pi);
	E = invertKTE(M,M,e);
	r_p = [a*(cos(E)-e);b*sin(E);0];
	v_p = (a*n/norm(r_p))*[-a*sin(E);b*cos(E);0];

	r = pCi'*r_p;
	v = pCi'*v_p;
end