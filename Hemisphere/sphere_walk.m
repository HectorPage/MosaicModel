function sphere_walk(epochs,initial_position,sigma, lower_scale_limit, upper_scale_limit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  GENERATING RANDOM WALK OVER SURFACE OF A HEMISPHERE  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Created by (and copyright held by) Hector Page on 14/10/2015

head_rotations = zeros(1,epochs+1);
head_rotations(1) = 0;

lower_scale_limit = abs(lower_scale_limit);   %set range of distance scaling factor - no <0 values
upper_scale_limit = abs(upper_scale_limit);


headings = zeros(epochs+1,3);
positions = zeros(epochs+1,3);




initial_heading = [0,0,0];
next_movement_direction = [0,0,0];

%%Reading in initial position
initial_position = initial_position/norm(initial_position); %normalised to keep on unit sphere
positions(1,:) = initial_position;  %recording the rat's start position

for counter = 1:epochs;
    
    %Randomly setting scaling factor for distance travelled.
    scaling = (upper_scale_limit-lower_scale_limit).*rand(1,1) + lower_scale_limit;
    
    %Setting initial position for this timestep
    initial_position = positions(counter,:);
    initial_position = initial_position/norm(initial_position); %normalised to keep on unit sphere
    
    %%Setting and recording initial heading - could randomise this?
    if initial_heading(:) == 0;
        %initial_heading = [initial_position(3),initial_position(3),...
         %   -(initial_position(1)) - initial_position(2)] * scaling;
        
         initial_heading = [0,1,0] * scaling;
         
        headings(1,:) = initial_heading;  
        movement_direction = initial_heading;
    else
        movement_direction = next_movement_direction;
    end
    
    
    movement_direction = movement_direction/norm(movement_direction);
    
    %%Calculating final position
    
    final_position = initial_position + (movement_direction * scaling);
    
    if final_position(3)< 0  %Avoiding rat walking off the edge of hemisphere  - maybe rethink?
        final_position(3) = 0.0005; %%If rat is at edge of hemisphere, C code doesn't work
    end
       
    final_position = final_position/norm(final_position); %normalised to keep on unit sphere.
    
    positions(counter+1,:) = final_position;
    
    %%Calculating unrotated_final_heading
    
    unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head
    
    
    %%Calculating next motion direction
    
    %To get the NEXT direction of movement: convert to quaternions to do rotation about the axis defined by
    %surface normal at final_position, and the angle defined by position_rotation
    
    if(initial_position(3) < 0.1 && headings(counter,3) < 0)  %Turn away from the edge of hemisphere - change params?
        
        if(initial_position(2)<0 && headings(counter,1)>0)
            position_rotation = normrnd(-45,sigma) * (pi/180);    %need to be in rad to do quat rotations
        end
        
        if(initial_position(2)<0 && headings(counter,1)<0)
            position_rotation = normrnd(45,sigma) * (pi/180); 
        end
        
        if(initial_position(2)>0 && headings(counter,1)<0)
            position_rotation = normrnd(-45,sigma) * (pi/180); 
        end
        
        if(initial_position(2)>0 && headings(counter,1)>0)
            position_rotation = normrnd(45,sigma) * (pi/180); 
        end
    else
        position_rotation = normrnd(0,sigma) * (pi/180);   %ditto, rad mode please
    end
    
    if position_rotation > pi
        position_rotation = pi - position_rotation;
    end
    
    head_rotations(counter+1) = position_rotation;
    
    %%Rotating heading
    surface_normal = final_position * 2; %This is the surface normal at initial position.
    
    next_movement_direction = rodrigues_rot(unrotated_final_heading,surface_normal,position_rotation);
   
    final_heading = next_movement_direction;
    
    headings(counter+1,:) = final_heading;
end




%%%NOW TO CORRECT THE HEADINGS TO LIE WITHIN THE TANGENT PLANE AT EACH
%%%POINT I.E. MAKING THEM ORTHOGONAL TO SURFACE NORMAL AT POSITION

corrected_headings = headings;
for index = 2:epochs  %first one is always ok, can't do last one, so ignore
 
    %This is a function I downloaded from the file exchange, it finds the
    %intersection of a plane defined by a point and a surface normal, with
    %a line defined by two points
    [intersection_point,~] = plane_line_intersect(2*positions(index,:),positions(index,:), [0 0 0],2*positions(index+1,:));

   
    %corrected heading is the vector from position(index) to the
    %intersection point
    
    corrected_headings(index,:) = intersection_point - positions(index,:);
end

headings = corrected_headings; %now making headings the corrected versions

%%%NOW PLOTTING HEADINGS%%%
figure('Name','Heading','NumberTitle','off');

hold on
quiver3(positions(:,1),positions(:,2),positions(:,3)...
    ,headings(:,1),headings(:,2),headings(:,3),'b');
hold on


view(-90,90);
axis([-1 1 -1 1 0 1]);

hold on

%Overlaying axes as vectors

axes_origin = zeros(6,3);
axes_vectors = [1,0,0;0,1,0;0,0,1;-1,0,0;0,-1,0;0,0,-1;]; %Providing origin of the axes

quiver3(axes_origin(:,1),axes_origin(:,2),axes_origin(:,3)...
    ,axes_vectors(:,1),axes_vectors(:,2),axes_vectors(:,3),0,'k');

hold on

plot3(positions(1,1),positions(1,2),positions(1,3),'X','MarkerSize',25,'MarkerFaceColor', 'r'); %Plotting start location

%Overlaying unit hemisphere to look nice/give context
sphere_radius = 1;
[sphere_x, sphere_y, sphere_z] = sphere(50);
sphere_x = sphere_x(26:end,:);       %# Keep top 26 x points
sphere_y = sphere_y(26:end,:);       %# Keep top 26 y points
sphere_z = sphere_z(26:end,:);       %# Keep top 26 z points
lightGrey = 0.9*[1 1 1]; % Making lines lighter, to not overwhelm vectors
surface(sphere_x,sphere_y,sphere_z,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
circle(0,0,sphere_radius);  %see code at end of function
axis square;
daspect([1 1 1]);

%%%NOW PLOTTING PATH%%%
figure('Name','Path','NumberTitle','off');
plot3(positions(:,1),positions(:,2),positions(:,3), 'r');
hold on

view(-90,90);
axis([-1 1 -1 1 0 1]);

hold on

%Overlaying axes as vectors
%Providing origin of the axes
axes_origin = zeros(6,3);
axes_vectors = [1,0,0;0,1,0;0,0,1;-1,0,0;0,-1,0;0,0,-1;];

quiver3(axes_origin(:,1),axes_origin(:,2),axes_origin(:,3)...
    ,axes_vectors(:,1),axes_vectors(:,2),axes_vectors(:,3),0,'k');

hold on

plot3(positions(1,1),positions(1,2),positions(1,3),'X','MarkerSize',25,'MarkerFaceColor', 'r');

%Overlaying unit sphere to look nice/give context
sphere_radius = 1;
[sphere_x, sphere_y, sphere_z] = sphere(50);
sphere_x = sphere_x(26:end,:);       %# Keep top 26 x points
sphere_y = sphere_y(26:end,:);       %# Keep top 26 y points
sphere_z = sphere_z(26:end,:);       %# Keep top 26 z points
lightGrey = 0.9*[1 1 1]; % Making lines lighter, to not overwhelm vectors
surface(sphere_x,sphere_y,sphere_z,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
circle(0,0,sphere_radius);
axis square;

daspect([1 1 1]);

%Writing path data for reading into c code.

extended_positions = zeros((epochs*100)+1,3);
extended_headings = zeros((epochs*100)+1,3);



for idx = 0:epochs-1
    for jdx = 1:100
        extended_positions((idx*100)+jdx,1) = positions(idx+1,1);
        extended_positions((idx*100)+jdx,2) = positions(idx+1,2);
        extended_positions((idx*100)+jdx,3) = positions(idx+1,3);
        
        extended_headings((idx*100)+jdx,1) = headings(idx+1,1);
        extended_headings((idx*100)+jdx,2) = headings(idx+1,2);
        extended_headings((idx*100)+jdx,3) = headings(idx+1,3);
    end
end

extended_positions((epochs*100)+1,1) = positions(epochs+1,1);
extended_positions((epochs*100)+1,2) = positions(epochs+1,2);
extended_positions((epochs*100)+1,3) = positions(epochs+1,3);

extended_headings((epochs*100)+1,1) = headings(epochs+1,1);
extended_headings((epochs*100)+1,2) = headings(epochs+1,2);
extended_headings((epochs*100)+1,3) = headings(epochs+1,3);

formatSpec = '%f\n';

fileID = fopen('positions_x.txt','w');
fprintf(fileID, formatSpec, extended_positions(:,1));
fclose(fileID);

fileID = fopen('positions_y.txt','w');
fprintf(fileID, formatSpec, extended_positions(:,2));
fclose(fileID);

fileID = fopen('positions_z.txt','w');
fprintf(fileID, formatSpec, extended_positions(:,3));
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen('headings_x.txt','w');
fprintf(fileID, formatSpec, extended_headings(:,1));
fclose(fileID);

fileID = fopen('headings_y.txt','w');
fprintf(fileID, formatSpec, extended_headings(:,2));
fclose(fileID);

fileID = fopen('headings_z.txt','w');
fprintf(fileID, formatSpec, extended_headings(:,3));
fclose(fileID);

fileID = fopen('head_rotations.txt','w');
fprintf(fileID,formatSpec, head_rotations*(180/pi));
fclose(fileID);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Keeping short form positions
fileID = fopen('positions_x_short.txt','w');
fprintf(fileID, formatSpec, positions(:,1));
fclose(fileID);

fileID = fopen('positions_y_short.txt','w');
fprintf(fileID, formatSpec, positions(:,2));
fclose(fileID);

fileID = fopen('positions_z_short.txt','w');
fprintf(fileID, formatSpec, positions(:,3));
fclose(fileID);

fileID = fopen('headings_x_short.txt','w');
fprintf(fileID, formatSpec, headings(:,1));
fclose(fileID);

fileID = fopen('headings_y_short.txt','w');
fprintf(fileID, formatSpec, headings(:,2));
fclose(fileID);

fileID = fopen('headings_z_short.txt','w');
fprintf(fileID, formatSpec, headings(:,3));
fclose(fileID);


%%%NOW VISUALISING A SMOOTHED VERSION OF THE PATH%%%

figure('Name','SMOOTH Path','NumberTitle','off')
fnplt(cscvn(positions.'),'r');
hold on
view(-90,90);
axis([-1 1 -1 1 0 1]);
hold on

%Overlaying axes as vectors
quiver3(axes_origin(:,1),axes_origin(:,2),axes_origin(:,3)...
    ,axes_vectors(:,1),axes_vectors(:,2),axes_vectors(:,3),0,'k');

hold on

plot3(positions(1,1),positions(1,2),positions(1,3),'X','MarkerSize',25,'MarkerFaceColor', 'r'); %Plotting start location

%Overlaying unit hemisphere to look nice/give context
sphere_radius = 1;
[sphere_x, sphere_y, sphere_z] = sphere(50);
sphere_x = sphere_x(26:end,:);       %# Keep top 26 x points
sphere_y = sphere_y(26:end,:);       %# Keep top 26 y points
sphere_z = sphere_z(26:end,:);       %# Keep top 26 z points
lightGrey = 0.9*[1 1 1]; % Making lines lighter, to not overwhelm vectors
surface(sphere_x,sphere_y,sphere_z,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
circle(0,0,sphere_radius);  %see code at end of function
axis square;
daspect([1 1 1]);

%%%NOW VISUALISING A SMOOTHED VERSION OF THE HEADINGS%%%

smoothed_path = fnplt(cscvn(positions.'));
smoothed_heading = fnplt(cscvn(headings.'));

figure('Name','SMOOTH Heading','NumberTitle','off')
quiver3(smoothed_path(1,:),smoothed_path(2,:),smoothed_path(3,:)...
    ,smoothed_heading(1,:),smoothed_heading(2,:),smoothed_heading(3,:));

hold on
grid off

view(-90,90);
axis([-1 1 -1 1 0 1]);

hold on

%Overlaying axes as vectors
quiver3(axes_origin(:,1),axes_origin(:,2),axes_origin(:,3)...
    ,axes_vectors(:,1),axes_vectors(:,2),axes_vectors(:,3),0,'k');

hold on

plot3(positions(1,1),positions(1,2),positions(1,3),'X','MarkerSize',25,'MarkerFaceColor', 'r'); %Plotting start location

%Overlaying unit hemisphere to look nice/give context
sphere_radius = 1;
[sphere_x, sphere_y, sphere_z] = sphere(50);
sphere_x = sphere_x(26:end,:);       %# Keep top 26 x points
sphere_y = sphere_y(26:end,:);       %# Keep top 26 y points
sphere_z = sphere_z(26:end,:);       %# Keep top 26 z points
lightGrey = 0.9*[1 1 1]; % Making lines lighter, to not overwhelm vectors
surface(sphere_x,sphere_y,sphere_z,'FaceColor', 'none','EdgeColor',lightGrey);
hold on
circle(0,0,sphere_radius);  %see code at end of function
axis square;

daspect([1 1 1]);

% figure()
% scatter3(smoothed_path(:,1),smoothed_path(:,2),smoothed_path(:,3));




%%PLOTTING HEATMAP%% - PROJECT TEMPORARILY ABANDONED, MIGHT RE-TRY

%converting (x,y,z) coords to sphericals (just phi and theta as rad is
%always 1)

%smoothed_circ_positions will be in format azimuth, elevation

% smoothed_circ_positions = zeros(epochs+1,2);
%
% smoothed_circ_positions(:,1) = atan2(smoothed_path(:,2),smoothed_path(:,1)) * (180/pi);    %azimuth
% smoothed_circ_positions(:,2) = atan2(smoothed_path(:,3),sqrt(smoothed_path(:,1).^2 + ...
%     smoothed_path(:,2).^2)) * (180/pi);    %elevation



% x = -180:12:180;                       % range for x bins
%
% if(x<0)
%     x = x+12/2;
% else
%     x = x-12/2;
% end
%
% y = x;  %y bins same as x bins
%
% figure()
% [N,C] = hist3(smoothed_circ_positions, 'edges', {x,y});
% %set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%
% [heat_x,heat_y,heat_z] = sphere(60);
% surface(heat_x,heat_y,heat_z,N);


end
function circle(x,y,r)  %this little function is to draw equator on at bottom of hemisphere
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%will be less smooth.
ang=0:0.01:2*pi;
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'k--');
end

function [I,check]=plane_line_intersect(n,V0,P0,P1)
%plane_line_intersect computes the intersection of a plane and a segment(or
%a straight line)
% Inputs: 
%       n: normal vector of the Plane 
%       V0: any point that belongs to the Plane 
%       P0: end point 1 of the segment P0P1
%       P1:  end point 2 of the segment P0P1
%
%Outputs:
%      I    is the point of interection 
%     Check is an indicator:
%      0 => disjoint (no intersection)
%      1 => the plane intersects P0P1 in the unique point I
%      2 => the segment lies in the plane
%      3=>the intersection lies outside the segment P0P1
%
% Example:
% Determine the intersection of following the plane x+y+z+3=0 with the segment P0P1:
% The plane is represented by the normal vector n=[1 1 1]
% and an arbitrary point that lies on the plane, ex: V0=[1 1 -5]
% The segment is represented by the following two points
% P0=[-5 1 -1]
%P1=[1 2 3]   
% [I,check]=plane_line_intersect([1 1 1],[1 1 -5],[-5 1 -1],[1 2 3]);

%This function is written by :
%                             Nassim Khaled
%                             Wayne State University
%                             Research Assistant and Phd candidate
%If you have any comments or face any problems, please feel free to leave
%your comments and i will try to reply to you as fast as possible.

I=[0 0 0];
u = P1-P0;
w = P0 - V0;
D = dot(n,u);
N = -dot(n,w);
check=0;
if abs(D) < 10^-7        % The segment is parallel to plane
        if N == 0           % The segment lies in plane
            check=2;
            return
        else
            check=0;       %no intersection
            return
        end
end

%compute the intersection parameter
sI = N / D;
I = P0+ sI.*u;

if (sI < 0 || sI > 1)
    check= 3;          %The intersection point  lies outside the segment, so there is no intersection
else
    check=1;
end
end