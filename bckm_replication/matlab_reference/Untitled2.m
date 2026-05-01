close all; clear all; clc;

load worktemp.mat

plot(worktemp.time(worktemp.bind:end),[worktemp.w.yt(worktemp.bind:end)*100 ...
    worktemp.w.mzy(worktemp.bind:end) worktemp.w.mly(worktemp.bind:end)...
    worktemp.w.mxy(worktemp.bind:end) worktemp.w.mgy(worktemp.bind:end)...
    ],'Linewidth',2);
legend('output','mz','ml','mx','mg','Location','best');
xlim([worktemp.time(worktemp.bind) worktemp.time(end)]);
title('Output and output components - US 1980.25:2014.5, 2008Q1 as base date');
print(1,'-djpeg','ycomponents.jpg');

plot(worktemp.time(worktemp.bind:end),[worktemp.w.ht(worktemp.bind:end)*100 ...
    worktemp.w.mzh(worktemp.bind:end) worktemp.w.mlh(worktemp.bind:end)...
    worktemp.w.mxh(worktemp.bind:end) worktemp.w.mgh(worktemp.bind:end)...
    ],'Linewidth',2);
legend('hours','mz','ml','mx','mg','Location','best');
xlim([worktemp.time(worktemp.bind) worktemp.time(end)]);
title('Hours and hours components - US 1980.25:2014.5, 2008Q1 as base date');
print(1,'-djpeg','hcomponents.jpg');

plot(worktemp.time(worktemp.bind:end),[worktemp.w.xt(worktemp.bind:end)*100 ...
    worktemp.w.mzx(worktemp.bind:end) worktemp.w.mlx(worktemp.bind:end) ...
    worktemp.w.mxx(worktemp.bind:end) worktemp.w.mgx(worktemp.bind:end) ...
    ],'Linewidth',2);
legend('investment','mz','ml','mx','mg','Location','best');
xlim([worktemp.time((worktemp.bind)) worktemp.time(end)]);
title('Investment and investment components - US 1980.25:2014.5');
print(1,'-djpeg','xcomponents.jpg');

plot(worktemp.time(worktemp.bind:end),[worktemp.w.zt(worktemp.bind:end)...
    worktemp.w.tault(worktemp.bind:end) worktemp.w.tauxt(worktemp.bind:end)...
    worktemp.w.gt(worktemp.bind:end)],'Linewidth',2);
legend('\omega_A','\omega_L','\omega_X','\omega_G','Location','best');
title('Wedges - US 1980.25:2014.5, 2008Q1 bse date');
xlim([worktemp.time(worktemp.bind) worktemp.time(end)]);
print(1,'-djpeg','wedges.jpg');

plot(worktemp.time(worktemp.bind:end),[worktemp.w.yt(worktemp.bind:end) ...
    worktemp.w.ht(worktemp.bind:end) worktemp.w.xt(worktemp.bind:end) ...
    worktemp.w.gt(worktemp.bind:end)],'Linewidth',2);
legend('output','hours','investment','government consumption','Location','best');
title('Observables - US 1980.25:2014.5');
xlim([worktemp.time(worktemp.bind) worktemp.time(end)]);
print(1,'-djpeg','observables.jpg');
