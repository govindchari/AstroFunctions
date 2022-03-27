function E = invertKTE(M,E0,e)
	E = E0;
	while (abs(M-E+e*sin(E)) > 1e-12)
		E = E - (M-E+e*sin(E))/(e*cos(E)-1);
	end
end