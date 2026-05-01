% clear all; close all; clc;

%% use this if wanting to run code for all the countries in the folder
dfiles = dir('d???.*');   %all data files 
ndcs   = size(dfiles,1);  %number of data files

for i = 1:ndcs
    dcnames{i} = dfiles(i).name(2:4);
end


%% data and estimation
mlestart = 1969.25;       %set starting period for MLE
mleend   = 2014.5;     %set ending period for MLE
freq     = 4;        %data frequency
adjc     = 2;        %1 for no, 2 for BGG, 3 for 4*BGG.


%% settings for plotting the wedges vs data
sy       = 1969;  %year
sq       =    1;  %quarter
ws       =  182;  %window size
wtitles  = {'\omega_A','\omega_L','\omega_X','\omega_G'};
hpcolor  = {'r','g','b','m'};
lw       = 1;     %linewidth for wedge plots


%% set uncmin parameters
nps = 10;    % number of uncmin runs
pb  = 0.99; % x0[k+1] = x0[k]*pb -> sequence of starting points for uncmin


% do bcas
for i = 1:ndcs
    load(['d',dcnames{i}]);                               
    worktemp.optimnum.nps = nps;
    worktemp.optimnum.pb  = pb;
    worktemp.mlestart     = mlestart;
    worktemp.mleend       = mleend;
    worktemp.time         = worktemp.mled(:,1);
    worktemp.adjc         = adjc;
    ind                   = find(worktemp.time==sy+sq/4);
    iniwd                 = worktemp.time(ind);
    worktemp.wsize        = [sy;sq;ws;ind;iniwd];
    worktemp.freq         = 4;
    save('worktemp.mat','worktemp','-mat');
    run('gmle.m'); 
    run('gwedges2.m');
    load worktemp.mat
    save(['bca',dcnames{i},'.mat'],'worktemp','-mat');
end


%% collect fi stats
bfiles = dir('bca???.*'); %all the bca exercises
nbcs   = size(bfiles,1);  %number of bca exercises
for i = 1:nbcs
    bcnames{i} = bfiles(i).name(4:6);
end

for i=1:nbcs
    load(['bca',bcnames{i}]);
    iniwd = worktemp.wsize(5);
    ind   = worktemp.wsize(4);
    wsize = worktemp.wsize(3);
    
    % errors for 1 wedge economies - output
    mzye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mzy(ind:ind+wsize-1)./worktemp.w.mzy(ind);
    mlye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mly(ind:ind+wsize-1)./worktemp.w.mly(ind);
    mxye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mxy(ind:ind+wsize-1)./worktemp.w.mxy(ind);
    mgye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mgy(ind:ind+wsize-1)./worktemp.w.mgy(ind);    
    worktemp.w.w1yerrors = [mzye mlye mxye mgye];
    temp1 = mean(worktemp.w.w1yerrors(:,1:4).^2,1);
    worktemp.w.w1yfz = (1/temp1(1))/sum(temp1.^-1);
    worktemp.w.w1yfl = (1/temp1(2))/sum(temp1.^-1);
    worktemp.w.w1yfx = (1/temp1(3))/sum(temp1.^-1);
    worktemp.w.w1yfg = (1/temp1(4))/sum(temp1.^-1);
    

    % errors for 1 wedge economies - hours
    mzhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mzh(ind:ind+wsize-1)./worktemp.w.mzh(ind);
    mlhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mlh(ind:ind+wsize-1)./worktemp.w.mlh(ind);
    mxhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mxh(ind:ind+wsize-1)./worktemp.w.mxh(ind);
    mghe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mgh(ind:ind+wsize-1)./worktemp.w.mgh(ind);    
    worktemp.w.w1herrors = [mzhe mlhe mxhe mghe];
    temp2 = mean(worktemp.w.w1herrors(:,1:4).^2,1);
    worktemp.w.w1hfz = (1/temp2(1))/sum(temp2.^-1);
    worktemp.w.w1hfl = (1/temp2(2))/sum(temp2.^-1);
    worktemp.w.w1hfx = (1/temp2(3))/sum(temp2.^-1);
    worktemp.w.w1hfg = (1/temp2(4))/sum(temp2.^-1);
    
    % errors for 1 wedge economies - investment
    mzxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mzx(ind:ind+wsize-1)./worktemp.w.mzx(ind);
    mlxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mlx(ind:ind+wsize-1)./worktemp.w.mlx(ind);
    mxxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mxx(ind:ind+wsize-1)./worktemp.w.mxx(ind);
    mgxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mgx(ind:ind+wsize-1)./worktemp.w.mgx(ind);    
    worktemp.w.w1xerrors = [mzxe mlxe mxxe mgxe];
    temp3 = mean(worktemp.w.w1xerrors(:,1:4).^2,1);
    worktemp.w.w1xfz = (1/temp3(1))/sum(temp3.^-1);
    worktemp.w.w1xfl = (1/temp3(2))/sum(temp3.^-1);
    worktemp.w.w1xfx = (1/temp3(3))/sum(temp3.^-1);
    worktemp.w.w1xfg = (1/temp3(4))/sum(temp3.^-1);
    
    % errors for 1 wedge economies - consumption
    mzce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mzc(ind:ind+wsize-1)./worktemp.w.mzc(ind);
    mlce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mlc(ind:ind+wsize-1)./worktemp.w.mlc(ind);
    mxce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mxc(ind:ind+wsize-1)./worktemp.w.mxc(ind);
    mgce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mgc(ind:ind+wsize-1)./worktemp.w.mgc(ind);    
    worktemp.w.w1cerrors = [mzce mlce mxce mgce];
    temp31 = mean(worktemp.w.w1cerrors(:,1:4).^2,1);
    worktemp.w.w1cfz = (1/temp31(1))/sum(temp31.^-1);
    worktemp.w.w1cfl = (1/temp31(2))/sum(temp31.^-1);
    worktemp.w.w1cfx = (1/temp31(3))/sum(temp31.^-1);
    worktemp.w.w1cfg = (1/temp31(4))/sum(temp31.^-1);
    
    % errors for 3 wedge economies - output
    mnozye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mnozy(ind:ind+wsize-1)./worktemp.w.mnozy(ind);
    mnolye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mnoly(ind:ind+wsize-1)./worktemp.w.mnoly(ind);
    mnoxye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mnoxy(ind:ind+wsize-1)./worktemp.w.mnoxy(ind);
    mnogye =   worktemp.w.Y(ind:ind+wsize-1,2)./worktemp.w.Y(ind,2) - ...
           worktemp.w.mnogy(ind:ind+wsize-1)./worktemp.w.mnogy(ind);    
    worktemp.w.w3yerrors = [mnozye mnolye mnoxye mnogye];
    temp4 = mean(worktemp.w.w3yerrors.^2,1);
    worktemp.w.w3yfz = (1/temp4(1))/sum(temp4.^-1);
    worktemp.w.w3yfl = (1/temp4(2))/sum(temp4.^-1);
    worktemp.w.w3yfx = (1/temp4(3))/sum(temp4.^-1);
    worktemp.w.w3yfg = (1/temp4(4))/sum(temp4.^-1);

    % errors for 3 wedge economies - hours
    mnozhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mnozh(ind:ind+wsize-1)./worktemp.w.mnozy(ind);
    mnolhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mnolh(ind:ind+wsize-1)./worktemp.w.mnoly(ind);
    mnoxhe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mnoxh(ind:ind+wsize-1)./worktemp.w.mnoxy(ind);
    mnoghe =   worktemp.w.Y(ind:ind+wsize-1,3)./worktemp.w.Y(ind,3) - ...
           worktemp.w.mnogh(ind:ind+wsize-1)./worktemp.w.mnogy(ind);    
    worktemp.w.w3herrors = [mnozhe mnolhe mnoxhe mnoghe];
    temp5 = mean(worktemp.w.w3herrors.^2,1);
    
    worktemp.w.w3hfz = (1/temp5(1))/sum(temp5.^-1);
    worktemp.w.w3hfl = (1/temp5(2))/sum(temp5.^-1);
    worktemp.w.w3hfx = (1/temp5(3))/sum(temp5.^-1);
    worktemp.w.w3hfg = (1/temp5(4))/sum(temp5.^-1);
    
    % errors for 3 wedge economies - investment
    mnozxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mnozx(ind:ind+wsize-1)./worktemp.w.mnozy(ind);
    mnolxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mnolx(ind:ind+wsize-1)./worktemp.w.mnoly(ind);
    mnoxxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mnoxx(ind:ind+wsize-1)./worktemp.w.mnoxy(ind);
    mnogxe =   worktemp.w.Y(ind:ind+wsize-1,4)./worktemp.w.Y(ind,4) - ...
           worktemp.w.mnogx(ind:ind+wsize-1)./worktemp.w.mnogy(ind);    
    worktemp.w.w3xerrors = [mnozxe mnolxe mnoxxe mnogxe];
    temp6 = mean(worktemp.w.w3xerrors.^2,1);
    
    worktemp.w.w3xfz = (1/temp6(1))/sum(temp6.^-1);
    worktemp.w.w3xfl = (1/temp6(2))/sum(temp6.^-1);
    worktemp.w.w3xfx = (1/temp6(3))/sum(temp6.^-1);
    worktemp.w.w3xfg = (1/temp6(4))/sum(temp6.^-1);
    
    % errors for 3 wedge economies - consumption
    mnozce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mnozc(ind:ind+wsize-1)./worktemp.w.mnozy(ind);
    mnolce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mnolc(ind:ind+wsize-1)./worktemp.w.mnoly(ind);
    mnoxce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mnoxc(ind:ind+wsize-1)./worktemp.w.mnoxy(ind);
    mnogce =   worktemp.w.Y(ind:ind+wsize-1,5)./worktemp.w.Y(ind,5) - ...
           worktemp.w.mnogc(ind:ind+wsize-1)./worktemp.w.mnogy(ind);    
    worktemp.w.w3cerrors = [mnozce mnolce mnoxce mnogce];
    temp7 = mean(worktemp.w.w3cerrors.^2,1);
    
    worktemp.w.w3cfz = (1/temp7(1))/sum(temp7.^-1);
    worktemp.w.w3cfl = (1/temp7(2))/sum(temp7.^-1);
    worktemp.w.w3cfx = (1/temp7(3))/sum(temp7.^-1);
    worktemp.w.w3cfg = (1/temp7(4))/sum(temp7.^-1);
    
    save(['bca',dcnames{i},'.mat'],'worktemp','-mat');  
end


    
    