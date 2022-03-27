function [v0,v] = lambert(r0,r,dt,tm)
    %Lambert solver using univeral variables and bisection

    %INPUTS
    %   r0 - initial heliocentric position 
    %   r - final heliocentric position
    %   dt - transfer time
    %   tm - binary variable (1 or -1) to indicate long way or short way
    %   transfer

    %OUTPUTS
    %   v0 - initial inertial velocity vector to complete transfer
    %   v - final inertial velocity vector

    mu = 2.9591220828559093e-4; %Sun's gravitational parameter
    flag = false;
	cn = dot(r0,r)/(norm(r0)*norm(r));
	sn = tm*sqrt(1-cn^2);
	A = tm*sqrt(norm(r)*norm(r0)*(1+cn));
	psin = 0.0;
	psinp1 = 0.0;
	psiup = 4.0*pi^2;
	psilow = -4.0*pi;
	chin = 0.0;
	c2 = (1/2);
	c3 = (1/6);
	yn = 0;
	dtn = 10000000000000.0;
	iter = 0;
	while (abs(dt-dtn) > 1e-6)
		iter = iter +1;
		yn = norm(r0) + norm(r) + (A*(psin*c3-1))/(sqrt(c2));
		if (A>0.0 && yn<0.0)
			psilow = 1.05*psin;
		else
			chin = sqrt(yn/c2);
			dtn = (chin^3*c3+A*sqrt(yn))/sqrt(mu);
			if (dtn <= dt)
				psilow = psin;
			else
				psiup = psin;
			end
			psinp1 = 0.5*(psiup+psilow);
			if (psinp1 >= 0)
				c2 = (1-cos(sqrt(psinp1)))/psinp1;
				c3 = (sqrt(psinp1)-sin(sqrt(psinp1)))/(psinp1^(1.5));
			else
				c2 = (1-cosh(sqrt(-psinp1)))/psinp1;
				c3 = (sinh(sqrt(-psinp1))-sqrt(-psinp1))/((-psinp1)^(1.5));
			end
		end
		if (iter>100)
            flag = true;
			break;
		end
		psin = psinp1;
	end
	f = 1-(yn/norm(r0));
	gd = 1-(yn/norm(r));
	g = A*sqrt(yn/mu);
	v0 = (r-f*r0)/g;
	v = (gd*r-r0)/g;
    if flag
        v0 = 100*ones(3,1);
        v = 100*ones(3,1);
    end
end
