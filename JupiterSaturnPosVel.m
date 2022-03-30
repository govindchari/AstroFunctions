function [r_e,r_v,v_e,v_v] = JupiterSaturnPosVel(t0,dt)

%INPUTS
%   t0 - Initial time (JD)
%   dt - Cruise duration (days)
%
%OUTPUTS
%   r_e - 3x1 col vector - Heliocentric position of Earth (in AU) at t0
%   r_v - 3x1 col vector - Heliocentric position of Venus (in AU) at t0+dt
%   v_e - 3x1 col vector - Heliocentric-Inertial velocity of Earth 
%                          (in AU/day) at t0
%   v_v - 3x1 col vector - Heliocentric-Inertial velocity of Venus 
%                          (in AU/day) at t0+dt

%Jupiter and Saturn orbital elements
aj = 5.202552135553255;
ej = 4.875572951666173e-2;
O_j = deg2rad(1.005161987950841e2);
I_j = deg2rad(1.303774910017153);
w_j = deg2rad(2.736952858336796e2);
tpj = 2459967.250901828986;
muj = 2.9619468438424372e-4;

as = 9.581741770422434;
es = 5.092477629223650e-2;
O_s = deg2rad(1.135931808557482e2);
I_s = deg2rad(2.489867712870254);
w_s = deg2rad(3.377632388714605e2);
tps = 2463572.349504575599;
mus = 2.9599678443364383e-4;


    [t_e,x_e_p,y_e_p,vx_e_p,vy_e_p,E_e] = kepler_2body(aj,ej,muj,tpj,t0);
    [t_v,x_v_p,y_v_p,vx_v_p,vy_v_p,E_v] = kepler_2body(as,es,mus,tps,t0);
    
    r_e_p=[x_e_p(1);y_e_p(1);0];
    v_e_p=[vx_e_p(1);vy_e_p(1);0];
    r_v_p=[x_v_p(length(vx_v_p));y_v_p(length(vx_v_p));0];
    v_v_p=[vx_v_p(length(vx_v_p));vy_v_p(length(vy_v_p));0];
    
    C1 = @(th) [1 0 0;0 cos(th) sin(th);0 -sin(th) cos(th)]; 
    C3 = @(th) [cos(th) sin(th) 0;-sin(th) cos(th) 0; 0 0 1];
    
    PCI_e=C3(w_j)*C1(I_j)*C3(O_j);
    PCI_v=C3(w_s)*C1(I_s)*C3(O_s);

    r_e =PCI_e'*r_e_p;
    r_v =PCI_v'*r_v_p;
    v_e =PCI_e'*v_e_p;
    v_v =PCI_v'*v_v_p;

function [t,x,y,vx,vy,E] = kepler_2body(a,e,mu,tp,ti) 
        % Two-body problem using Newton-Raphson inversion of Kepler's
        % time equation
        % Inputs:
        %   a - semi-major axis
        %   e - eccentricity
        %   mu - gravitational parameter (with same distance units as a)
        %
        % Output:
        %   t - 1000 element time array between 0 and the orbital period
        %   x,y - Orbital radius components for one orbit in the perifocal
        %         frame (e,q) directions. Same size as t.
        %   E - Specific orbital energy over one orbital period. Same
        %       size as t.

        l=a*(1-e^2);
        n = sqrt(mu/a^3);       %mean motion
        h = sqrt(mu*l);       %angular momentum
        Tp = 2*pi*sqrt(a^3/mu);     %orbital period

        %create time array from 0 to To with 1000 total points
        t = linspace(ti,ti+dt,1e3).';

        %calculate the Mean anomaly at each of the values in t
        M = n*(t-tp); %mean anomaly
        M=mod(M,2*pi);

        %Newton-Raphson to find eccentric anomaly. Iterate to 
        %machine-percision tolerance tol=eps(2*pi)
        tol = eps(2*pi);

        for i=1:length(M)
            if (M(i)/(1-e))<sqrt((6*(1-e))/e)
                E0=M(i)/(1-e);
            else
                E0=(6*M/e)^(1/3);
            end
            E1=E0;
            E2=E1-(M(i)-E1+e*sin(E1))/(e*cos(E1)-1);
            while (abs(E2-E1)>tol)
                E1=E2;
                E2=E1-(M(i)-E1+e*sin(E1))/(e*cos(E1)-1);
            end
            Ecc(i)=E1;
        end

        Eccd=n./(1-e*cos(Ecc));

        %calculate x and y positions in the perifocal frame using Kepler's
        %equations
        b= a*sqrt(1-e^2);
        x = a*(cos(Ecc)-e)     ;
        y = b*sin(Ecc)     ;

        %calcualte the x and y velocities in the perifocal frame
        vx = -a*Eccd.*sin(Ecc)    ;
        vy = b*Eccd.*cos(Ecc)    ;

        %calculate the specific orbital energy
        E = 0.5*(vx.^2+vy.^2) - mu./sqrt(x.^2+y.^2) ;

        %force all outputs to be col vectors
        x = x(:);
        y = y(:);
        E = E(:);

    end


%ensure outputs are column vectors
r_e = r_e(:);
r_v = r_v(:);
v_e = v_e(:);
v_v = v_v(:);
end