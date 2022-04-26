function [t,z,E] = n_body_integration(mu_list, z0, tspan, dt)
    %Energy scaled by a factor of G
    num_bodies = length(z0)/6;
    t = [0:dt:tspan];
    z = zeros(length(z0),length(t));
    E = zeros(length(t),1);
    z(:,1) = z0;
    E(1) = n_body_energy(z(:,1));
    for k=2:length(t)
        z(:,k) = propagate_leapfrog(@n_body_force,z(1:3*num_bodies,k-1),z(3*num_bodies+1:6*num_bodies,k-1),dt);
        E(k) = n_body_energy(z(:,k));
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