%Change working directory to where VIIRS data are stored
%
%/projectnb/modislc/users/emelaas/VIIRS/NPP_VDNE/
data_dir = 'R:/Mann Research/India Night Time Lights/May2015/';
d1 = dir('NPP_VDNE*hdf');
d2 = dir('NPP_CMIP*hdf');
d3 = dir('NPP_VMAE*hdf');
d1


%Loop through each DNB/Cloud Mask/Surface Reflectance file
for i = 1:length(d1)
    disp(i)
    cd(data_dir)

    
    %Set path names
    path1 = [data_dir d1(i).name];
    path2 = [data_dir d2(i).name];
    path3 = [data_dir d3(i).name];
    
    %Read in data from HDFs
    Radiance = hdfread(path1,'VIIRS_EV_DNB_SDR', 'Fields', 'Radiance');
    Latitude = hdfread(path1,'VIIRS_EV_DNB_SDR', 'Fields', 'Latitude');
    Longitude = hdfread(path1,'VIIRS_EV_DNB_SDR', 'Fields', 'Longitude');
    QF3_VIIRSCMIP = hdfread(path2,'VCM_IP', 'Fields', 'QF3_VIIRSCMIP');
    Latitude2 = hdfread(path3,'VIIRS_EV_750M_SDR', 'Fields', 'Latitude');
    Longitude2 = hdfread(path3,'VIIRS_EV_750M_SDR', 'Fields', 'Longitude');
    
    %Create new HDF5 files
    h5create([path1(64:76) '.h5'],'/Radiance',size(Radiance));
    h5create([path1(64:76) '.h5'],'/Latitude',size(Latitude));
    h5create([path1(64:76) '.h5'],'/Longitude',size(Longitude));
    h5create([path1(64:76) '.h5'],'/CloudMask',size(QF3_VIIRSCMIP));
    h5create([path1(64:76) '.h5'],'/Latitude2',size(Latitude2));
    h5create([path1(64:76) '.h5'],'/Longitude2',size(Longitude2));
    
    %Write data to HDF5 file
    h5write([path1(64:76) '.h5'], '/Radiance', Radiance);
    h5write([path1(64:76) '.h5'], '/Latitude', Latitude);
    h5write([path1(64:76) '.h5'], '/Longitude', Longitude);
    h5write([path1(64:76) '.h5'], '/CloudMask', QF3_VIIRSCMIP);
    h5write([path1(64:76) '.h5'], '/Latitude2', Latitude2);
    h5write([path1(64:76) '.h5'], '/Longitude2', Longitude2);
end
