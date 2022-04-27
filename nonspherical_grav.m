function f = nonspherical_grav(r,Rb,mu,l_list, m_list, Clm_bar_list, Slm_bar_list)
    % r is position vector in ECEF Frame 
    % Rb is body equitorial radius
    % r, mu, and Rb must have consistent units
    % f is specific perturbing force in ECEF

    phi = acos(r(3)/norm(r));
    th = atan2(r(2),r(1));
    [Clm_list, Slm_list] = unnormalize_legendre(l_list, m_list, Clm_bar_list, Slm_bar_list);

    l_start = l_list(1);
    l_end = l_list(length(l_list)-1);

    dR_dr = 0;
    dR_dth = 0;
    dR_dphi = 0;

    for l=l_start:l_end
        Plm_list = legendre(l,cos(phi));
        Plmom_list = legendre(l-1,cos(phi));
        rBrl = (Rb/norm(r))^l;
        for m=0:l
            idx = (l^2+l-4)/2 + m;
            Plm = Plm_list(m+1);
            if m<l
                Plmom = Plmom_list(m+1);
            else
                Plmom = 0;
            end
            Clm = Clm_list(idx);
            Slm = Slm_list(idx);
            dR_dr = dR_dr + rBrl*(l+1)*Plm*(Clm*cos(m*th)+Slm*(sin(m*th)));
            dR_dth = dR_dth + rBrl*m*Plm*(-Clm*sin(m*th)+Slm*cos(m*th));
            dR_dphi = dR_dphi + rBrl*(l*cos(phi)*Plm-(l+m)*Plmom)*((Clm*cos(m*th)+Slm*sin(m*th))/(sin(phi)));
        end
    end

    dR_dr = -(mu/norm(r)^2) * dR_dr;
    dR_dth = (mu/norm(r))*dR_dth;
    dR_dphi = (mu/norm(r))*dR_dphi;

    rho = norm(r(1:2));
    f1 = (r(1)/norm(r))*(dR_dr+(r(3)/(norm(r)*rho))*dR_dphi) - (r(2)/rho^2)*dR_dth;
    f2 = (r(1)/rho^2)*dR_dth + (r(2)/norm(r))*(dR_dr + (r(3)/(norm(r)*rho))*dR_dphi);
    f3 = (r(3)/norm(r))*dR_dr - (rho/norm(r)^2)*dR_dphi;
    f = [f1;f2;f3];
end