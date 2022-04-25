function znext = propagate_euler(f, t, z, dt)
    znext = z + f(t, z) * dt;
end