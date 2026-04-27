function dd = runnit()


% countries
dcnames = {'AUS','AUT','BEL','CAN','DNK','FIN','FRA','GER',...
           'IRL','ICE','ISR','ITA','JPN','KOR','LUX','MEX','NLD',...
           'NOR','NZL','SPA','SWE','SWI','UK','USA'};
       
for i = 1:size(dcnames,2)
    cd ..;
    cd(dcnames{i});
    dirsource = cd;
    cd ..;
    cd('AAA-Perfect Foresight adj');
    dirtarget = cd;
    copyfile(dirsource,dirtarget);
    run('bca_steady2.m');
    run('bca_wedges2.m');
    run('bca_simul2.m');
    copyfile(dirtarget,dirsource);

end
    