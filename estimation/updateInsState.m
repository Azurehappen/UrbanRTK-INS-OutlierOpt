function x_plus = updateInsState(p, dx, x_prior, Rot_e2g)

if p.state_mode == p.ins_mode
    x_plus = x_prior;
    % correct pos and vel states
    x_plus(1:6) = x_prior(1:6) + Rot_e2g(1:6,1:6)' * dx(1:6);
    % Correct the quaternion (Farrell Appendix.D Exercise D.4)
    rho = dx(7:9);
    rho_bar = 0.5*rho;
    % q_u= sqrt(1-(norm(rho_bar))^2);
    % q_e=[q_u;rho_bar];
    if norm(rho_bar) > 1
        q_e = 1/sqrt(1+norm(rho_bar)^2)*[1; rho_bar];
    else
        q_e = [sqrt(1-norm(rho_bar)^2); rho_bar];
    end
    x_plus(7:10) = quatInv(quatMult(q_e,quatInv(x_prior(7:10,1))));
    % Correct other states
    x_plus(11:end) = x_prior(11:end) + dx(10:length(x_prior)-1);

else
    x_plus = x_prior + Rot_e2g(1:9,1:9)'*dx(1:9);
end

end