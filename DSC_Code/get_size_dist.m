%  Get_traj() reads the fortran binary file containing the trajectory
%  profiles for each scenario

traj_td = zeros(size(spec_list,2),6,ndays);

for dd = 1:ndays
for spec = spec_list
    ii = strmatch(spec, spec_list,'exact');
    
    %READ DATA
    fname = ['D:/Trajectory/' cell2mat(scenario.age) '/' ...
        cell2mat(scenario.tag) '/' cell2mat(scenario.loc) '/' cell2mat(spec)...
        '.' cell2mat(scenario.loc) '.' cell2mat(scenario.date(dd)) '.trj'];
    [day,hour,conc1] = ...
        textread(fname, '%f %f %*s %*s %f %*f %*f %*f %*f %*f %*f %*f %*f %*f');
    
    weed = isnan(conc1);
    conc1(weed) = 0;
    traj_td(ii,1:6,dd) = conc1;

end
end
