load modsims; load wedges;

figure
plot((2008.25:0.25:2011.75)',[wedges(:,1)./wedges(1,1) wedges(:,2)./wedges(1,2)]);
title('Perfect foresight vs standard BCA - Efficiency wedge 2008Q1-2011Q3');
legend('Perfect foresight','standard BCA');
print(1,'PFvsBCA - Awedge','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[wedges(:,3)./wedges(1,3) wedges(:,4)./wedges(1,4)]);
title('Perfect foresight vs standard BCA - Labor wedge 2008Q1-2011Q3');
legend('Perfect foresight','standard BCA');
print(1,'PFvsBCA - Lwedge','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[wedges(:,5)./wedges(1,5) wedges(:,6)./wedges(1,6)]);
title('Perfect foresight vs standard BCA - Investment wedge 2008Q1-2011Q3');
legend('Perfect foresight','standard BCA');
print(1,'PFvsBCA - Xwedge','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[wedges(:,7)./wedges(1,7) wedges(:,8)./wedges(1,8)]);
title('Perfect foresight vs standard BCA - government wedge 2008Q1-2011Q3');
legend('Perfect foresight','standard BCA');
print(1,'PFvsBCA - Gwedge','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[datvars(:,1)./datvars(1,1) modsims(:,1)./modsims(1,1) modsims(:,2)./modsims(1,2)]);
title('Perfect foresight vs standard BCA - A_t models 2008Q1-2011Q3');
legend('Output','Perfect foresight','standard BCA');
print(1,'PFvsBCA - Amodel','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[datvars(:,1)./datvars(1,1) modsims(:,3)./modsims(1,3) modsims(:,4)./modsims(1,4)]);
title('Perfect foresight vs standard BCA - \tau_L models 2008Q1-2011Q3');
legend('Output','Perfect foresight','standard BCA');
print(1,'PFvsBCA - Lmodel','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[datvars(:,1)./datvars(1,1) modsims(:,5)./modsims(1,5) modsims(:,6)./modsims(1,6)]);
title('Perfect foresight vs standard BCA - \tau_x models 2008Q1-2011Q3');
legend('Output','Perfect foresight','standard BCA');
print(1,'PFvsBCA - Xmodel','-dpng');close 1;

figure
plot((2008.25:0.25:2011.75)',[datvars(:,1)./datvars(1,1) modsims(:,7)./modsims(1,7) modsims(:,8)./modsims(1,8)]);
title('Perfect foresight vs standard BCA - gt models 2008Q1-2011Q3');
legend('Output','Perfect foresight','standard BCA');
print(1,'PFvsBCA - Gmodel','-dpng');close 1;