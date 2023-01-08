function out = hep_npm(t,in)

% Model of Ca2+ dynamics during hormone stimulation
% using an open cell model with SOCE

c = in(1);
h = in(2);
p = in(3);

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       A K_PLC ...
       gamma ct;
   
% Functions

ce = gamma.*(ct - c);

PLC = A.*c.^2./(K_PLC^2 + c.^2);

phi_c = c.^4./(c.^4 + Kc^4);
phi_p = p.^2./(p.^2 + Kp^2);
phi_p_down = Kp^2./(p.^2 + Kp^2);

h_inf = Kh^4./(Kh^4 + c.^4);
tau = tau_max.*K_tau^4./(K_tau^4 + c.^4);

Jserca = Vs.*(c.^2 - Kbar.*ce.^2)./(c.^2 + Ks^2);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + Kb.*(beta + alpha));

Jipr = Kf.*Po.*(ce-c);

% DE's of the model

out(1) = Jipr - Jserca;
out(2) = (h_inf - h)./tau;
out(3) = tau_p.*(PLC - p);
out = out';

end