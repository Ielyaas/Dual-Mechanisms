function out = hep_fixed_h(~,in,param_ct)

% Function file for fixed h

c = in(1);
p = in(2);

ce = param_ct.gamma.*(ct - c);

PLC = param_ct.R_act.*c.^2./(param_ct.K_PLC^2 + c.^2);

phi_c = c.^4./(c.^4 + param_ct.Kc^4);
phi_p = p.^2./(p.^2 + param_ct.Kp^2);
phi_p_down = param_ct.Kp^2./(p.^2 + param_ct.Kp^2);
h_inf = param_ct.Kh^4./(param_ct.Kh^4 + c.^4);

Jserca = param_ct.Vs.*(c.^2 - param_ct.Kbar.*ce.^2)./(c.^2 + param_ct.Ks^2);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + param_ct.Kb.*(beta + alpha));

Jipr = param_ct.Kf.*Po.*(ce-c);

% DE's of the model

out(1) = Jipr - Jserca;
out(2) = param_ct.tau_p*(PLC - p);

end