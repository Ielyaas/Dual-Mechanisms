function out = hep_open(t,in)

% Model of Ca2+ dynamics during hormone stimulation in an open cell

c = in(1);
ce = in(2);
h = in(3);
p = in(4);

global Kf Kb Kc Kh Kp ...
       tau_max K_tau ...
       Vs Kbar Ks ...
       tau_p ...
       V_PLC K_PLC ...
       Vpm Kpm ...
       alpha0 alpha1 ...
       delta gamma Vp_start;
   
% Functions

if t < Vp_start
    
    delta = 2.5;
else
    delta = 0;
    
end

m_a = c.^4./(Kc^4 + c.^4);
m_b = m_a;
h_a = Kh^4./(Kh^4 + c.^4);
h_inf = h_a;
b = p.^2./(Kp^2 + p.^2);
a = 1- b;
Beta = b.*m_b.*h;
Alpha = a.*(1-(m_a.*h_a));

Po = Beta./(Beta + Kb.*(Beta + Alpha));

tau_h = tau_max.*(K_tau^4./(K_tau^4 + c.^4));

J_ipr = Kf.*Po.*(ce - c);

Jserca = Vs.*((c.^2 - Kbar.*ce.^2)./(Ks^2 + c.^2));

nu_PLC = V_PLC.*c.^2./(K_PLC^2 + c.^2);

Jpm = Vpm.*c.^2./(Kpm^2 + c.^2);
Jin = alpha0 + alpha1*V_PLC;

% DE's of the model

out(1) = J_ipr - Jserca + delta.*(Jin - Jpm);
out(2) = gamma.*(Jserca - J_ipr);
out(3) = (h_inf - h)./tau_h;
out(4) = tau_p.*(nu_PLC - p);
out = out';

end