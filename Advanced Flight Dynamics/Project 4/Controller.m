function U = Controller(u)
% Evaluate the sample's LQ tracking law using the current reference input.
x = u(1:6);
yd = u(7:12);

[R, B, P, g_aux] = Riccati(yd);
U = R \ (B' * (g_aux - P*x));
end



