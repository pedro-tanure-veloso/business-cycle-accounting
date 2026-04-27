mzy = worktemp.w.mzy(ind:ind+wsize-1)./worktemp.w.mzy(ind);
mly =  worktemp.w.mly(ind:ind+wsize-1)./worktemp.w.mly(ind);
mxy = worktemp.w.mxy(ind:ind+wsize-1)./worktemp.w.mxy(ind);
mgy = worktemp.w.mgy(ind:ind+wsize-1)./worktemp.w.mgy(ind);
y = worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2);

plot([(worktemp.w.mzy+worktemp.w.mly+worktemp.w.mxy+worktemp.w.mgy)/4]) ...
    worktemp.w.Y(:,2)])


plot([Y(:,1)-Y(Y0,1) YMz(:,1)-YMz(Y0,1)+...
                     YMl(:,1)-YMl(Y0,1)+...
                     YMx(:,1)-YMx(Y0,1)+...
                     YMg(:,1)-YMg(Y0,1)])

                 
y = worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2);
X = [worktemp.w.mzy(ind:ind+wsize-1)./worktemp.w.mzy(ind) ...
     worktemp.w.mly(ind:ind+wsize-1)./worktemp.w.mly(ind) ...
     worktemp.w.mxy(ind:ind+wsize-1)./worktemp.w.mxy(ind) ...
     worktemp.w.mgy(ind:ind+wsize-1)./worktemp.w.mgy(ind)];
 
n = size(X,1); 
[b,bint,r,rint,stats] = regress(y-mean(y),X-repmat(mean(X,1),[n 1]));

[zb,zbint,zr,zrint,zstats] = regress(y,[ones(n,1) X(:,1)]);
[lb,lbint,lr,lrint,lstats] = regress(y,[ones(n,1) X(:,2)]);
[xb,xbint,xr,xrint,xstats] = regress(y,[ones(n,1) X(:,3)]);
[gb,gbint,gr,grint,gstats] = regress(y,[ones(n,1) X(:,4)]);

fzi = zstats(1)/(zstats(1)+lstats(1)+xstats(1)+gstats(1));
fli = lstats(1)/(zstats(1)+lstats(1)+xstats(1)+gstats(1));
fxi = xstats(1)/(zstats(1)+lstats(1)+xstats(1)+gstats(1));
fgi = gstats(1)/(zstats(1)+lstats(1)+xstats(1)+gstats(1));

[b,resnorm] = lsqnonneg(X(:,1),y)