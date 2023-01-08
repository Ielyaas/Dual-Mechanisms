clear all
close all
clc

% Open probability of IPR

% Parameters
% IPR2
Kh = 0.08; Kc = 0.5; Kp = 0.1; Kb = 0.3;

% IPR1
Kh1 = 0.11; Kc1 = 0.9; Kp1 = 0.3; Kb1 = 0.5;


p = 10;
c = 0:0.01:5;

% Functions
% IPR2
m_a = c.^4./(Kc^4 + c.^4);
m_b = m_a;
h_a = Kh^4./(Kh^4 + c.^4);
h_inf = h_a;
h = h_inf;
b = p.^2./(Kp^2 + p.^2);
a = 1- b;
Beta = b.*m_b.*h;
Alpha = a.*(1-(m_a.*h_a));

Po = Beta./(Beta + Kb.*(Beta + Alpha));

% IPR1
m_a1 = c.^4./(Kc1^4 + c.^4);
m_b1 = m_a1;
h_a1 = Kh1^4./(Kh1^4 + c.^4);
h_inf1 = h_a1;
h1 = h_inf1;
b1 = p.^2./(Kp1^2 + p.^2);
a1 = 1- b1;
Beta1 = b1.*m_b1.*h1;
Alpha1 = a1.*(1-(m_a1.*h_a1));

Po1 = Beta1./(Beta1 + Kb1.*(Beta1 + Alpha1));

% Plots

figure(1)
plot(c,Po,'LineWidth',1);
hold on
plot(c,Po1,'LineWidth',1);
legend({'IPR2','IPR1'}, 'Fontsize',18);
xlabel('Ca^{2+} \muM','Fontsize',20,'Fontweight','b')
ylabel('Open probability, P_o','Fontsize',20,'Fontweight','b')