function out = hep_J(t,in)

% Model of Ca2+ dynamics during hormone stimulation
% using an open cell model

c = in(1);
h = in(2);
ce = in(3);
p = in(4);

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       A K_PLC ...
       Vpm Kpm ...
       tau_i alpha0 alpha1 ...
       delta gamma pm_block in_block;
   
% Functions

%ce = gamma.*(ct - c);

if t < pm_block
    
    delta = 2.5;
else
    delta = 0;
    
end

if t < in_block
    
    tau_i = 1;
else
    tau_i = 0;
    
end

PLC = A.*c.^2./(K_PLC^2 + c.^2);

phi_c = c.^4./(c.^4 + Kc^4);
phi_p = p.^2./(p.^2 + Kp^2);
phi_p_down = Kp^2./(p.^2 + Kp^2);

Jpm = Vpm.*c.^2./(Kpm^2 + c.^2);
Jin = tau_i*(alpha0 + alpha1*A);

h_inf = Kh^4./(Kh^4 + c.^4);
tau = tau_max*K_tau^4./(K_tau^4 + c.^4);

Jserca = Vs.*(c.^2 - Kbar.*ce.^2)./(c.^2 + Ks^2);

beta = phi_p.*phi_c.*h;
alpha = phi_p_down.*(1-phi_c.*h_inf);

Po = beta./(beta + Kb.*(beta + alpha));

Jipr = Kf.*Po.*(ce-c);

% DE's of the model

out(1) = Jipr - Jserca + delta.*(Jin - Jpm);
out(2) = (h_inf - h)./tau;
out(3) = gamma.*(Jserca - Jipr);
out(4) = tau_p.*(PLC - p);
out = out';

end