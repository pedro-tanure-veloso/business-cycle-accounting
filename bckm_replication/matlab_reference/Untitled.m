plot([worktemp.w.yt(100:140)*100 worktemp.w.mzy(100:140) worktemp.w.mly(100:140) worktemp.w.mgy(100:140)],'LineWidth',2); 

legend('output','\omega_A','\omega_L','\omega_G');

title('CD 2008b - Obs 100 to 200');xlim([1 100]);

w1t = [worktemp.w.mzy worktemp.w.mly worktemp.w.mxy worktemp.w.mgy];