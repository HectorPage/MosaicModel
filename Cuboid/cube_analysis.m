function  cube_analysis(path_time,tsize)

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

Pvectors = zeros(1,numel(InputPvectors)+1); %This stores the input to the network (i.e. combined increment)
Pvectors = HD_Pvectors;
%Pvectors(2:end) = InputPvectors;
Pvectors = Pvectors * (pi/180);


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

head_azimuth = zeros(1,path_timesteps);     %azimuth of current heading rotated up to earth horiztonal
difference_measure = zeros(1,path_timesteps);


exp_Ncell_dir = zeros(1,path_timesteps);


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
correction_factor = atan2(sin(start_facing_azimuth-Initial_Pvector), cos(start_facing_azimuth-Initial_Pvector));

calculated_final_vector = zeros(path_timesteps,3);
rot_product = zeros(1,path_timesteps);

wall_heading_angle = zeros(1,path_timesteps); %hplds heading in current wall reference frame (rel to up or north)

TopIndex = false(1,path_timesteps);
NorthIndex = false(1,path_timesteps);
SouthIndex = false(1,path_timesteps);
EastIndex = false(1,path_timesteps);
WestIndex = false(1,path_timesteps); %These are logical indexes to display only data from each wall.

cell_azimuth = zeros(1,path_timesteps);
  
for idx = 1:path_timesteps %only run for path timesteps, not path timesteps +1, as last element is not used in simulation

    %Setting surfnorm based on current cuboid face
    if (pos_z(idx) == 1.6)
        surfnorm = [0 0 1];
        
        TopIndex(idx) = 1;
        
    elseif(pos_y(idx) == -0.5) %South wall
        surfnorm = [0 -1 0];
        
        SouthIndex(idx) = 1;
        
    elseif(pos_y(idx) == 0.5)   %North wall
        surfnorm = [0 1 0];
      
        NorthIndex(idx) = 1;
        
    elseif(pos_x(idx) == -0.5) %West wall
        surfnorm = [-1 0 0];
      
        WestIndex(idx) = 1;
        
    elseif(pos_x(idx) == 0.5)  %East wall
        surfnorm = [1 0 0];
        
        EastIndex(idx) = 1;
        
    end
    
    [rot_axis, rot_angle] = euler_rotation(surfnorm,[0 0 1]);   %find how to align surfnorm with z axis.
    
    head = [head_x(idx) head_y(idx) head_z(idx)];
    head = head/(norm(head));
    
    
    
    if(all(surfnorm==[0 0 1]))
        calculated_final_vector(idx,:) = head;
    else
        calculated_final_vector(idx,:) = rodrigues_rot(head,rot_axis,-rot_angle);   %do the rotation
        calculated_final_vector(idx,:) = calculated_final_vector(idx,:)/(norm(calculated_final_vector(idx,:)));
    end
    
    rot_product(idx) = dot(calculated_final_vector(idx,:),[0 0 1]);
    head_azimuth(idx) = atan2(calculated_final_vector(idx,2),calculated_final_vector(idx,1));
    
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
    
    
     exp_Ncell_dir(idx) = (cell_current_azimuth - head_azimuth(idx)) * (180/pi);
    
    if(exp_Ncell_dir(idx)<0)
        exp_Ncell_dir(idx) = exp_Ncell_dir(idx) + 360.0;
    end
    
    
    distance1 = abs(cell_current_azimuth - head_azimuth(idx));
    distance2 = (2*pi)-distance1;
    
    if distance1<distance2
        difference_measure(idx) = distance1;
    else
        difference_measure(idx) = distance2;
    end
    
    difference_measure(idx) = difference_measure(idx) * (180/pi);
 

    if(head_z(idx) == 0)
        %wall_heading_angle(idx) = atan2(head(2),head(1)); %take azimuth angle of heading
        wall_heading_angle(idx) = acos(dot([-1 0 0],head));
        
        if(dot(surfnorm,cross(head,[-1 0 0])) < 0) %signing angle
             wall_heading_angle(idx) = -wall_heading_angle(idx);
    
         end
        
        
    else
        wall_heading_angle(idx) = acos(dot([0 0 1],head));%angle from UP
        
         if(dot(surfnorm,cross(head,[0 0 1])) < 0) %signing angle
             wall_heading_angle(idx) = -wall_heading_angle(idx);
    
         end
    end
    
    if(wall_heading_angle(idx)<0)
        wall_heading_angle(idx) = wall_heading_angle(idx) + (2*pi); %putting it in the correct range
    end
    
    cell_azimuth(idx) = cell_current_azimuth;
    
end


%Loading in 2D positions
xid = fopen('pos2d_x.txt','r');
pos2d_x = fscanf(xid, '%f');
fclose(xid);
clear xid;

yid = fopen('pos2d_y.txt','r');
pos2d_y = fscanf(yid, '%f');
fclose(yid);
clear yid;

    new_NorthIndex = false(1,path_timesteps/100);
    new_TopIndex = false(1,path_timesteps/100);
    new_SouthIndex = false(1,path_timesteps/100);
    new_EastIndex = false(1,path_timesteps/100);
    new_WestIndex = false(1,path_timesteps/100);
    
    new_difference_measure = zeros(1,path_timesteps/100);
    new_exp_Ncell_dir = zeros(1,path_timesteps/100);
    new_wall_heading_angle = zeros(1,path_timesteps/100);
    
    new_cell_azimuth = zeros(1,path_timesteps/100);
    
for idx = 1:(path_timesteps/100)
    new_difference_measure(idx) = difference_measure((idx*100));    
    new_exp_Ncell_dir(idx) = exp_Ncell_dir((idx*100));    
    new_wall_heading_angle(idx) = wall_heading_angle((idx*100));
    
    new_NorthIndex(idx) = logical(NorthIndex((idx*100)));
    new_TopIndex(idx) = logical(TopIndex((idx*100)));
    new_SouthIndex(idx) = logical(SouthIndex((idx*100)));
    new_EastIndex(idx) = logical(EastIndex((idx*100)));
    new_WestIndex(idx) = logical(WestIndex((idx*100)));
    
    new_cell_azimuth(idx) = cell_azimuth((idx*100));
end




%%%Plotting here the preferred firing direction of a single cell on each
%%%wall

%Load firing rates

%HDRates = zeros(500, path_timesteps);
HDRates = file_load(500, path_timesteps, 'HDRates.bdat');


x = 1:path_timesteps/100;
y1 = 1:500;
[X1, Y1] = meshgrid(x, y1);

figure()
h1 = surf(X1, Y1, HDRates(:,(1:path_timesteps/100)));



short_rates = zeros(500,path_timesteps/100);

difference = zeros(1,path_timesteps/100);


for idx = 1:(path_timesteps/100)
    short_rates(:,idx) = HDRates(:,(idx*100));
    
    new_cell_azimuth(idx) = pi-(new_cell_azimuth(idx));
    
    if(new_cell_azimuth(idx)<0)
        new_cell_azimuth(idx) = new_cell_azimuth(idx) + (2*pi);
    end
    
    difference1 = abs(new_cell_azimuth(idx)-new_wall_heading_angle(idx));
    difference2 = (2*pi) - difference1;
    
    if(difference1<difference2)
        difference(idx) = difference1;
    else
        difference(idx) = difference2;
    end
end




figure()
plot(new_wall_heading_angle * (180/pi));
hold on
plot((new_cell_azimuth * (180/pi)),'r');

figure();
plot(difference*(180/pi));

figure()
polar(new_wall_heading_angle(new_TopIndex),HDRates(250,new_TopIndex));

%this measure above works

figure()
scatter(pos2d_x(1:path_timesteps/100),pos2d_y(1:path_timesteps/100),45,difference*(180/pi),'Filled');
colormap(jet);
colorbar;
caxis([1 180])
hold on
line([-2 -2],[0 2.6],'color', 'r', 'LineStyle', '--');
line([-2 -1],[2.6 2.6],'color', 'r', 'LineStyle', '--');
line([-2 2],[1.6 1.6], 'color','r', 'LineStyle', '--');
line([-1 -1],[0 2.6], 'color','r', 'LineStyle', '--');
line([0 0],[0 1.6],'color', 'r', 'LineStyle', '--');
line([1 1],[0 1.6], 'color','r', 'LineStyle', '--');
title('Shift In Preferred Firing Direction');

% figure()
% scatter(pos2d_x(1:path_timesteps/100),pos2d_y(1:path_timesteps/100),45,new_wall_heading_angle,'Filled');
% colormap(jet);
% colorbar;
% caxis([1 360])
% hold on
% line([-2 -2],[0 2.6],'color', 'r', 'LineStyle', '--');
% line([-2 -1],[2.6 2.6],'color', 'r', 'LineStyle', '--');
% line([-2 2],[1.6 1.6], 'color','r', 'LineStyle', '--');
% line([-1 -1],[0 2.6], 'color','r', 'LineStyle', '--');
% line([0 0],[0 1.6],'color', 'r', 'LineStyle', '--');
% line([1 1],[0 1.6], 'color','r', 'LineStyle', '--');

figure();
plot(new_difference_measure);
xlim([0,path_timesteps/100]);
title('Error in path integration');
%ylim([0 180]);





figure();
scatter(pos2d_x(1:path_timesteps/100),pos2d_y(1:path_timesteps/100),45,new_difference_measure,'Filled');
colormap(jet);
colorbar;
caxis([0 180])
title('Error in path integration');
hold on
line([-2 -2],[0 2.6],'color', 'r', 'LineStyle', '--');
line([-2 -1],[2.6 2.6],'color', 'r', 'LineStyle', '--');
line([-2 2],[1.6 1.6], 'color','r', 'LineStyle', '--');
line([-1 -1],[0 2.6], 'color','r', 'LineStyle', '--');
line([0 0],[0 1.6],'color', 'r', 'LineStyle', '--');
line([1 1],[0 1.6], 'color','r', 'LineStyle', '--');




%%Now doing Kate's analysis

new_spacing = [360.0
                      337.5
                      315
                      292.5
                      270
                      247.5
                      225
                      202.5
                      180
                      157.5
                      135
                      112.5
                      90
                      67.5
                      45
                      22.5
                      0];
           
       
colours = [218, 37, 128
          176, 10, 128
          127, 0, 128
           79, 10, 128
           37, 37 128
           10, 79, 128
           0, 128, 128
           10, 176, 128
           37, 218, 128
           79, 245, 128
           128, 255, 128
           176, 245, 128
           218, 218, 128
           245, 176, 128
           255, 128, 128
           245, 79, 128
           218, 37, 128]./255;
           
         
       
       katescolourmap = interp1(new_spacing/360,colours,linspace(0,1,255));
%        
%        I = linspace(0,1,360);
%        imagesc(I(ones(1,1),:)');
%        colormap(katescolourmap);
%        
     
       
    
figure();
 annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Direction Indicated By North Cell', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center','FontSize', 20);
subplot(1,2,1);
scatter(pos2d_x(1:6000),pos2d_y(1:6000),45,new_exp_Ncell_dir,'Filled');
view([0 90]);
grid off;
set(gca,'XTickLabel','','YTickLabel','')
set(gca,'TickLength',[ 0 0 ])
daspect([1,1,1]);
caxis([0 360]);
hold on
line([-2 -2],[0 2.6],'color', 'k', 'LineStyle', '--');
line([-2 -1],[2.6 2.6],'color', 'k', 'LineStyle', '--');
line([-2 2],[1.6 1.6], 'color','k', 'LineStyle', '--');
line([-1 -1],[0 2.6], 'color','k', 'LineStyle', '--');
line([0 0],[0 1.6],'color', 'k', 'LineStyle', '--');
line([1 1],[0 1.6], 'color','k', 'LineStyle', '--');
set(gca,'visible','off')


direction_names = {'NW', 'W', 'SW', 'S', 'SE', 'E', 'NE', 'N'};
       index = 0;
       
       x=ones(1,360);
       subplot(1,2,2);
       h=pie(x);
       colormap(katescolourmap);
    for hc=1:2:length(h); % segments
        set(h(hc),'EdgeColor','none');
    end
    for hc=2:2:length(h); % texts
        if mod(hc,45)~=0;
            delete(h(hc)); % delete some texts
        else
            set(h(hc),'string',direction_names(index+1), 'FontSize', 10); % set new texts
            index = index+1;
        end
    end      
    
    g = get(gca, 'pos');
    g(1) = g(1) - 0.08;
    g(2) = g(2) + 0.58;
    g(3) = g(3)/3;
    g(4) = g(4)/3;
    set(gca,'pos',g);



end

function rates = file_load(cells, steps, fname)
        
       rates = zeros(cells, steps);
       
       fid = fopen(fname, 'rb');
       
       rates = fread(fid, [steps, cells], 'float32')';
       
       fclose(fid);
       
  
end
