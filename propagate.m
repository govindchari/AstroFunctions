function [r,v] = propagate(r0,v0,t,mu)
	[a,e,Omega,I,omega,tp] = rv2kepler(r0,v0,mu);
	CO = [cos(Omega) sin(Omega) 0;-sin(Omega) cos(Omega) 0;0 0 1];
	CI = [1 0 0;0 cos(I) sin(I);0 -sin(I) cos(I)];
	Co = [cos(omega) sin(omega) 0;-sin(omega) cos(omega) 0;0 0 1];
	pCi = Co*CI*CO;
	
	b = a*sqrt(1-e^2);
	n = sqrt(mu/a^3);
	M = mod(n*(t-tp),2*pi);

	E = invertKTE(M,M,e);
	r_p = [a*(cos(E)-e);b*sin(E);0];
	v_p = (a*n/norm(r_p))*[-a*sin(E);b*cos(E);0];

	r = pCi'*r_p;
	v = pCi'*v_p;

end