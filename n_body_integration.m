function [tret,z,E] = n_body_integration(mu_list, z0, tspan, dt, trec)
    %Energy scaled by a factor of G
    %trec is the recording interval must be greater or equal to dt
    num_bodies = length(z0)/6;
    t = [0:dt:tspan];
    z = zeros(length(z0),round(length(t)/trec));
    E = zeros(round(length(t)/trec),1);
    z(:,1) = z0;
    ztemp = z0;
    E(1) = n_body_energy(z(:,1));
    for k=2:length(t)
        ztemp = propagate_leapfrog(@n_body_force,ztemp(1:3*num_bodies),ztemp(3*num_bodies+1:6*num_bodies),dt);
        Etemp = n_body_energy(ztemp);
        if (abs(round(t(k)/trec)*trec-t(k)) <= 0.5*dt)
            idx = round(t(k)/trec);
            z(:,idx) = ztemp;
            E(idx) = Etemp;
            tret(idx) = t(k);
        end
        if (mod(k,round(length(t)/10)) == 0)
            sprintf('Percent Complete: %g', 10*(k/round(length(t)/10)))
        end
    end
    function f = n_body_force(q)
        f = zeros(3*num_bodies,1);
        for j=1:num_bodies
            force = zeros(3,1);
            for i=1:num_bodies
            if (i~=j)
                mu = mu_list(i);
                ri = q(3*i-2:3*i);
                rj = q(3*j-2:3*j);
                rij = ri-rj;
                force = force + (mu/norm(rij)^3)*rij;
            end
            end
            f(3*j-2:3*j) = force;
        end
    end
    function E = n_body_energy(z)
        q = z(1:3*num_bodies);
        qd = z(3*num_bodies+1:6*num_bodies);
        T = 0;
        U = 0;
        for j=1:num_bodies
            T = T + 0.5 * mu_list(j) * dot(qd(3*j-2:3*j),qd(3*j-2:3*j));
            for i=j+1:num_bodies
                ri = q(3*i-2:3*i);
                rj = q(3*j-2:3*j);
                rij = ri-rj;
                U = U - (mu_list(i) * mu_list(j))/norm(rij);
            end
        end
        E = (T+U);
    end
end