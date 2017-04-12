function  kate_analysis(path_time,tsize)

%Path time in seconds

path_timesteps = uint32(path_time/tsize);

%Loading Pvectors%
fid = fopen('Pvector.bdat','rb');
HD_Pvectors = fread(fid, path_timesteps, 'float32');
fclose(fid);
clear fid;

%Loading in Input PVector
fid = fopen('InputPvector.bdat', 'rb');
InputPvectors = fread(fid, path_timesteps,'float32');
fclose(fid);
clear fid;

Pvectors = zeros(1,numel(InputPvectors)+1);
Pvectors(1) = HD_Pvectors(1);
Pvectors = InputPvectors;
Pvectors = Pvectors * (pi/180);


%Loading in increments
fid = fopen('IdealIncrement.bdat', 'rb');
IdInc = fread(fid, path_timesteps,'float32');
fclose(fid);
clear fid;

fid = fopen('ActualIncrement.bdat', 'rb');
AcInc = fread(fid, path_timesteps,'float32');
fclose(fid);
clear fid;

%Loading in headings
xid = fopen('headings_x.txt','rb');
head_x = fscanf(xid, '%f');
fclose(xid);
clear xid;

yid = fopen('headings_y.txt','r');
head_y = fscanf(yid, '%f');
fclose(yid);
clear yid;

zid = fopen('headings_z.txt','r');
head_z = fscanf(zid, '%f');
fclose(zid);
clear zid;


%Loading in positions
xid = fopen('positions_x.txt','r');
pos_x = fscanf(xid, '%f');
fclose(xid);
clear xid;

yid = fopen('positions_y.txt','r');
pos_y = fscanf(yid, '%f');
fclose(yid);
clear yid;

zid = fopen('positions_z.txt','r');
pos_z = fscanf(zid, '%f');
fclose(zid);
clear zid;

%Loading in reference vectors
fid = fopen('refvecs.bdat','rb');
refvecs = fread(fid, [path_timesteps, 3] ,'double');
fclose(fid);
clear fid;

%Loading in ROTATED reference vectors
fid = fopen('rotated refvecs.bdat','rb');
ROTrefvecs = fread(fid, [path_timesteps, 3] ,'double');
fclose(fid);
clear fid;

%Loading in ROTATED headings
fid = fopen('rotated headings.bdat','rb');
ROTheadings = fread(fid, [path_timesteps, 3] ,'double');
fclose(fid);
clear fid;

%Loading start angles%
fid = fopen('startangle.bdat','rb');
startangle = fread(fid, path_timesteps, 'double');
fclose(fid);
clear fid;


head_azimuth = zeros(1,path_timesteps);     %azimuth of current heading projected to earth horiztonal
difference_measure = zeros(1,path_timesteps);


%Calculating what azimuth rat is facing when HD activity is in original position, using that to
%get the Pvector that represents North

start_facing_azimuth = atan2(head_y(1),head_x(1));
if start_facing_azimuth < 0
    start_facing_azimuth = start_facing_azimuth + (2*pi); %putting in 0 to 2pi range
end

%Need to find correction factor to convert Pvector to azimuth:
%1. Find start facing azimuth
%2. Find distance to current azimuth from Pvector

Initial_Pvector = HD_Pvectors(1) * (pi/180);
%correction_factor = atan2(sin(start_facing_azimuth-Pvectors(1)), cos(start_facing_azimuth-Pvectors(1)));
correction_factor = atan2(sin(start_facing_azimuth-Initial_Pvector), cos(start_facing_azimuth-Initial_Pvector));


for idx = 1:path_timesteps %only run for path timesteps, not path timesteps +1, as last element is not used in simulation

    head_azimuth(idx) = atan2(head_y(idx),head_x(idx));
    
    if head_azimuth(idx) < 0
        head_azimuth(idx) = head_azimuth(idx) + (2*pi);
    end
    
    cell_current_azimuth = Pvectors(idx) + correction_factor;
    
    if cell_current_azimuth >(2*pi)
        cell_current_azimuth = cell_current_azimuth - (2*pi);
    end
    
    if cell_current_azimuth <0 
        cell_current_azimuth = cell_current_azimuth + (2*pi);
    end
    
    distance1 = abs(cell_current_azimuth - head_azimuth(idx));
    distance2 = (2*pi)-distance1;
    
    if distance1<distance2
        difference_measure(idx) = distance1;
    else
        difference_measure(idx) = distance2;
    end
    
    difference_measure(idx) = difference_measure(idx) * (180/pi);
 
    
end


figure();
plot(difference_measure);
xlim([0,6000]);
title('Error in path integration');
%ylim([0 180]);

figure();
scatter3(pos_x(1:6000),pos_y(1:6000),pos_z(1:6000),45,difference_measure,'Filled');
colormap(jet);
colorbar;
caxis([0 180])
title('Error in path integration');


new_difference_measure = difference_measure;
new_difference_measure(new_difference_measure<40) = 40;
new_difference_measure(new_difference_measure>40) = 180;

AbsIncDev = zeros(1,6000);
AbsIncDevMod = zeros(1,6000);
AbsIncDev(1:6000) = abs(abs(IdInc)-abs(AcInc));
AbsIncDevMod(AbsIncDev<40) = 40;
AbsIncDevMod(AbsIncDev>40) = 160;

% figure();
% scatter3(pos_x(3910:3920),pos_y(3910:3920),pos_z(3910:3920),45,AbsIncDevMod(3910:3920),'Filled');
% colormap(1-hot);
% colorbar;
% title('Deviation of increment from ideal');
% caxis([0 180]);

figure();
scatter3(pos_x(1:6000),pos_y(1:6000),pos_z(1:6000),45,new_difference_measure,'Filled');
colormap(1-hot);
colorbar;
caxis([0 180])
title('High error locations');

figure();
scatter3(pos_x(1:6000),pos_y(1:6000),pos_z(1:6000),45,AbsIncDevMod,'Filled');
colormap(jet);
colorbar;
caxis([0 180])
title('High increment errors');

%%Checking Increment%%
figure()
plot(IdInc,'b');
hold on
plot(AcInc,'g');
legend('Ideal Increment','Actual Increment');

figure()
plot(abs(abs(IdInc)-abs(AcInc)));
title('Deviation of increment from ideal');

%%Kate's slope of plane check%%
inclination = zeros(1,path_timesteps);
RVdot = zeros(1,path_timesteps);

for idx = 1:path_timesteps
    surface_normal = [pos_x(idx)*2 pos_y(idx)*2 pos_z(idx)*2];
    surface_normal = surface_normal/(norm(surface_normal));
    
    RV = [refvecs(idx,1) refvecs(idx,2) refvecs(idx,3)];
    RV = RV/(norm(RV));
    
    RVdot(idx) = dot(surface_normal,RV);
    
    inclination(idx) = acos(dot(surface_normal,[0 0 1])) * (180/pi); 
    
    if(startangle(idx)>180)
        startangle(idx) = startangle(idx) - 360.0;
    end
    
end

increment_error = abs(abs(IdInc)-abs(AcInc));

figure()
scatter(inclination, increment_error)
xlabel('Tangent Angle');
ylabel('Increment Error');

compoundCondInd = (inclination >= 80) & (inclination <= 90);

figure()
scatter(abs(startangle(compoundCondInd)), increment_error(compoundCondInd));
xlabel('Starting Angle Deviation');
ylabel('Increment Error');
title('Tangent inclined between 80 and 90 degrees');




end
