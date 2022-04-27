function [Clm_arr, Slm_arr] = unnormalize_legendre(l_arr, m_arr, Clm_bar_arr, Slm_bar_arr)
    %Unnormalizes Clm and Slm coefficients for nonspherial gravitational
    %potential
    len = length(l_arr);
    Clm_arr = zeros(len,1);
    Slm_arr = zeros(len,1);
    for i=1:len
        l = l_arr(i);
        m = m_arr(i);
        num = factorial(l-m)*(2*l+1)*(2-eq(m,0));
        den = factorial(l+m);
        Clm_arr(i) = sqrt((num)/(den))*Clm_bar_arr(i);
        Slm_arr(i) = sqrt((num)/(den))*Slm_bar_arr(i);
    end
end