% Define your transfer function G(s) = num(s)/den(s)
num = [-173 264]; 
den = [1 24 -400];          % example: s^2 + 4 s + 8
G = tf(num,den);

[GM, PM, Wcg, Wcp] = margin(G);     % GM (abs), PM (deg), Wcg/Wcp (rad/s)
margin(G)                           % plots Bode with margins marked

% If there are multiple crossings (e.g., nonminimum-phase), use:
am = allmargin(G);                  % returns ALL gain/phase margins