function H = MagResp(z,z1,z2,p1,K)
    z1c = conj(z1);
    z2c = conj(z2);
    p1c = conj(p1);
    H = K*(z-z1).*(z-z1c).*(z-z2).*(z-z2c)./(z-p1)./(z-p1c)./(z.^2);
end
