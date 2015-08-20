%Change working directory to where VIIRS data are stored
data_dir = '/groups/manngroup/India VIIRS/2015/';
cd(data_dir)
d1 = dir('NPP_VDNE*hdf');

%Loop through each DNB/Cloud Mask/Surface Reflectance file
moon_phase_angle = NaN*zeros(length(d1),1);
moon_illum_frac = NaN*zeros(length(d1),1);
year = NaN*zeros(length(d1),1);
doy = NaN*zeros(length(d1),1);
time = NaN*zeros(length(d1),1);
for i = 1:length(d1)
    disp(i)
    
    %Set path name
    path1 = [data_dir d1(i).name];
    
    %Extract date and time
    year(i,1) = str2double(path1(1,49:52));
    doy(i,1) = str2double(path1(1,53:55));
    time(i,1) = str2double(path1(1,57:60));
    
    %Read in Moon data from HDFs
    t = hdfinfo(path1);
    moon_phase_angle(i,1) = t.Attributes(4).Value;
    moon_illum_frac(i,1) = t.Attributes(5).Value;
end

%Concatenate data columns and write new matrix to CSV
data = [year doy time moon_illum_frac moon_phase_angle];
csvwrite('moon_info.csv',data);
