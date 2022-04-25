function [t,z] = n_body_integration(mu_list, z0, tspan, dt)
    num_bodies = length(z0)/6;
    t = [0:dt:tspan];
    z = zeros(length(z0),length(t));
    z(:,1) = z0;
    for k=2:length(t)
        z(:,k) = propagate_leapfrog(@solar_system_force,z(1:3*num_bodies,k-1),z(3*num_bodies+1:6*num_bodies,k-1),dt);
    end
    function f = solar_system_force(q)
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
end