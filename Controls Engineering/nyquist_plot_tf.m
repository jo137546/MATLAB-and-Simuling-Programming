s = tf('s');
G = tf([-173 264],[1 24 -400]);

Ks = [0.1 0.5 1 2 5];                 % try a few gains
opts = nyquistoptions; 
opts.ShowFullContour = 'on';          % show ±jω halves
figure; hold on
for k = Ks
    nyquistplot(k*G, opts);           % Nyquist of L(s)=K*G(s)
end
plot(-1,0,'rx','MarkerSize',10,'LineWidth',2)  % the -1 point
axis equal; grid on
legend(arrayfun(@(k)sprintf('K = %.2g',k),Ks,'uni',0), 'Location','best');
title('Nyquist of L(s)=K G(s)   (neg. unity feedback)');

opts = nyquistoptions;
opts.ShowFullContour = 'on';          % include the mirror (−jω) half
figure; nyquistplot(G, {1e-3, 1e3}, opts); hold on
plot(-1,0,'rx','MarkerSize',10,'LineWidth',2)  % critical point
axis equal; grid on; title('Nyquist of G(j\omega)');