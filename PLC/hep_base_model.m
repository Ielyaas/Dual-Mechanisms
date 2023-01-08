function out = hep_base_model(t,in,param)

% Base hepatocyte model

c = in(1);
ct = in(2);
h = in(3);
p = in(4);

% Functions

% if t < 20
%     L = 0;
% elseif t < 150 
%     L = 0.05;
% elseif t < 250
%     L = 0.1;
% else
%     L = param.V_plc*c.^2./(c.^2 + param.Kplc^2);
% end

% if t < 250
%     param.tau_p = 0;
% else
%     param.tau_p = 1;
% end
% 
% if t < 250
%     param.Vplc = 0;
% else 
%     param.Vplc = 0.4;
% end

ce = param.gamma.*(ct - c);

phi_c = c.^4./(c.^4 + param.Kc^4);
phi_p = p.^2./(p.^2 + param.Kp^2);
phi_p_down = param.Kp^2./(p.^2 + param.Kp^2);

h_inf = param.Kh^4./(param.Kh^4 + c.^4);
tau = param.tau_max.*param.K_tau^4./(param.K_tau^4 + c.^4);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + param.Kb.*(beta + alpha));

L = param.V_plc*c.^2./(c.^2 + param.Kplc^2);

Jipr = param.Kf.*Po.*(ce - c);
Jserca = param.Vs.*(c.^2 - param.Kbar.*ce.^2)./(c.^2 + param.Ks^2);

Jin = param.alpha0 + param.alpha1*param.Kce^4./(ce.^4 + param.Kce^4);
Jpm = param.Vpm.*c.^2./(c.^2 + param.Kpm^2);

% DE's of the model

out(1) = (Jipr - Jserca + param.delta.*(Jin - Jpm))/param.tau_c;
out(2) = param.delta.*(Jin - Jpm);
out(3) = (h_inf - h)./tau;
out(4) = param.tau_p.*(L - p);
out = out';

end

