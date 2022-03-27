function [a,e,Omega,I,omega,tp] = rv2kepler(r,v,mu) 

%Convert orbital position and velocity to orbital elements
%Inputs:
%   r - 3x1 array: orbital position vector
%   v - 3x1 array: orbital velocity vector
%   mu - scalar: gravitational parameter
%
%Outputs:
%   a - scalar: semi-major axis (infinite for parabolic orbits)
%   e - scalar: eccentricity
%   I - scalar: inclination (in radians, range 0 to pi)
%   omega - scalar: argument of periapsis (in radians, range 0 to 2*pi)
%   Omega - scalar: longitude of the ascending node (in radians, range 0 to
%                   2*pi)
%   tp - scalar: time of periapse passage, assuming that time of 
%                observation is zero. For closed orbits, this should
%                be a strictly positive value.

%All input units are assumed to be self-consistent.

%ensure that the inputs are column vectors
    r = r(:);
    v = v(:);
    
    A=eye(3);
    e1=A(:,1);
    e2=A(:,2);
    e3=A(:,3);
    
    %Calculates intermediate vectors
    h=cross(r,v);
    e_v=(cross(v,h))/mu-(r/norm(r));
    nhat=cross(e3,h)/norm(cross(e3,h));  
    
    %Calculates Geometric Parameters
    a = (dot(h,h))/(mu*(1-dot(e_v,e_v))); %semi-major axis
    e = norm(e_v); %eccentricity
    I = acos(h(3)/norm(h)); %inclination
    omegasign=dot(cross(nhat,e_v),h)/norm(dot(cross(nhat,e_v),h));
    omega = atan2(omegasign*(norm(cross(nhat,e_v))/(norm(e_v))),(dot(nhat,e_v)/(norm(e_v)))); %argument of periapsis
    Omegasign=dot(cross(e1,nhat),e3)/norm(dot(cross(e1,nhat),e3));
    Omega = atan2(Omegasign*norm(cross(e1,nhat)), nhat(1) ); %longitude of ascending node
    
    %Calculates time values (E,tp)
    tol=1e-6;
    nusign=dot(cross(e_v,r),h)/norm(dot(cross(e_v,r),h));
    nu= atan2(nusign*norm(cross(r,e_v)/((norm(r)*norm(e_v)))),dot(e_v,r)/((norm(r)*norm(e_v))));
    
    if abs(e-1)<tol
        a=1/0;
        E=tan(nu/2);
        l=(2*norm(r))/(1+E^2);
        np=2*sqrt(mu/l^3);
        tp=-(E+E^3/3)/np;
    elseif e<1
        nu=mod(nu,2*pi);
        E=atan2(norm(r)*sin(nu)/(a*sqrt(1-e^2)),(a*e+norm(r)*cos(nu))/(a));
        n=sqrt(mu/a^3);
        tp=(e*sin(E)-E)/n;
        period=2*pi*sqrt(a^3/mu);
        if tp<0
            tp=period+tp;
        end
    elseif e>1
        E=asinh(-(norm(r)*sin(nu))/(a*sqrt(e^2-1)));
        nh=sqrt(-mu/a^3);
        tp=(E+e*((r*sin(nu))/(a*sqrt(e^2-1)))/nh);
    end
    %ensure that all relevant angle outputs are in 0, 2*pi
    omega = mod(omega,2*pi);
    Omega = mod(Omega,2*pi);
end