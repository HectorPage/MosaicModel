function cuboid_walk(epochs,initial_position,sigma, lower_scale_limit, upper_scale_limit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     GENERATING RANDOM WALK OVER SURFACE OF A CUBOID    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Created by (and copyright held by) Hector Page on 23/11/2015




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%The script initially creates the path and headings on a the cuboid%%
%%flattened out onto a 2D surface, then constructs this in 3D after.%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


head_rotations = zeros(1,epochs+1);
head_rotations(1) = 0;

lower_scale_limit = abs(lower_scale_limit);   %set range of distance scaling factor - no <0 values
upper_scale_limit = abs(upper_scale_limit);


headings = zeros(epochs+1,2);
positions = zeros(epochs+1,2);



initial_heading = [0,0];
next_movement_direction = [0,0];

%%Reading in initial position

positions(1,:) = initial_position;  %recording the rat's start position

for counter = 1:epochs;
    
    %Randomly setting scaling factor for distance travelled.
    scaling = (upper_scale_limit-lower_scale_limit).*rand(1,1) + lower_scale_limit;
    
    %Setting initial position for this timestep
    initial_position = positions(counter,:);
    
    
    %%Setting and recording initial heading - could randomise this?
    if initial_heading(:) == 0;       
         initial_heading = ([-1,0]/(norm([-1 0]))) * scaling;
        headings(1,:) = initial_heading;  
        movement_direction = initial_heading;
    else
        movement_direction = next_movement_direction;
    end

    movement_direction = movement_direction/(norm(movement_direction));
    
    %%Calculating final position
    
    final_position = initial_position + (movement_direction * scaling);
    
    unrotated_final_heading = (final_position - initial_position); %final heading if rat had not rotated head    
    
    
    %%Now stopping rat from moving beyond cuboid edges%%

    %%Transitioning from the top to other walls
     if final_position(2) > 1.6 && initial_position(1) <-1  
         
       if final_position(2) > 2.6 %going over top down east wall - FIXED
         thestring = ['Top to E: Timestep',num2str(counter)];
           disp(thestring);
          final_position(1) = 0 + abs(final_position(1)) - 1;           
          final_position(2) = (1.6 - (final_position(2)-2.6));
          unrotated_final_heading = [-movement_direction(1) -movement_direction(2)];
       
       
     elseif final_position(1) <-2 %going over top down north wall - FIXED
            thestring = ['Top to N: Timestep',num2str(counter)];
           disp(thestring);
           fp1 = final_position(1); 
           fp2 = final_position(2); 
           final_position(1) = 1 + (2.6-fp2);
           final_position(2) = 1.6 - (abs(fp1) - 2);
           unrotated_final_heading = [movement_direction(1) -movement_direction(2)];
      
       
       elseif final_position(1) > -1 %going over top down south wall - FIXED
           thestring = ['Top to S: Timestep',num2str(counter)];
           disp(thestring);
           fp1 = final_position(1);
           fp2 = final_position(2);
           final_position(1) = -1 + (fp2-1.6);
           final_position(2) = 1.6 - (fp1 +1);
           unrotated_final_heading = [unrotated_final_heading(2) -unrotated_final_heading(1)];
       end
     end
    
     %%Transitioning between North and West walls
    if final_position(2) < 0  %At bottom, don't go underneath (cuboid on a surface...)
         thestring = ['Avoid floor: Timestep',num2str(counter)];
           disp(thestring);
        final_position(2) = 0; 
    end
   
    if final_position(1) < -2 %FIXED West to north
         thestring = ['W to N: Timestep',num2str(counter)];
           disp(thestring);
        final_position(1) = final_position(1)+4;
    end
    
    if final_position(1)>2 %FIXED North to west
         thestring = ['N to W: Timestep',num2str(counter)];
           disp(thestring);
        final_position(1) = final_position(1)-4;
    end
       
    
    %%Getting from other walls to the top
    if final_position(2) > 1.6 && final_position(1) >-1 %%getting to top surface
        if(initial_position(1) < 0)
            %south to top
            thestring = ['S to top: Timestep',num2str(counter)];
            disp(thestring);
            fp1 = final_position(1);
            fp2 = final_position(2);
            final_position(1) = -1 -(fp2-1.6);
            final_position(2) = 1.6 + fp1+1;
            unrotated_final_heading = [-unrotated_final_heading(2) unrotated_final_heading(1)]; %90 deg rot
        elseif(initial_position(1) < 1)
            %east to top - FIXED
             thestring = ['E to top: Timestep',num2str(counter)];
           disp(thestring);
            final_position(1) = -1 - final_position(1);
            final_position(2) = 2.6 - (final_position(2) - 1.6);
            unrotated_final_heading = -unrotated_final_heading; %180 deg rot
        else
            %north to top - FIXED
             thestring = ['N to top: Timestep',num2str(counter)];
           disp(thestring);
            fp1 = final_position(1);
            fp2 = final_position(2);
            final_position(1) = -2 + (fp2 - 1.6);
            final_position(2) = 2.6 - (fp1-1);
            unrotated_final_heading = [unrotated_final_heading(2) -unrotated_final_heading(1)]; %90 deg rot
        end
       
    end
    

    positions(counter+1,:) = final_position;
    
    %%Calculating unrotated_final_heading
    
    
    unrotated_final_heading = unrotated_final_heading/(norm(unrotated_final_heading)) * scaling;
     
    %%Calculating next motion direction
    
    %To get the NEXT direction of movement:  do rotation
    if(initial_position(2)< 0.1 && headings(counter,2)<0)
        if(headings(counter,1)<0)
            position_rotation = normrnd(-45,sigma) * (pi/180); 
        else
            position_rotation = normrnd(45,sigma) * (pi/180);
        end
    else
     position_rotation = normrnd(0,sigma) * (pi/180);   %ditto, rad mode please
    end
    
    if position_rotation > pi
        position_rotation = pi - position_rotation;
    end
    
    head_rotations(counter+1) = position_rotation;

    final_heading = [((cos(position_rotation)*unrotated_final_heading(1))-(sin(position_rotation)*unrotated_final_heading(2)))...
              , ((sin(position_rotation)*unrotated_final_heading(1))+(cos(position_rotation)*unrotated_final_heading(2)))];
   
          final_heading = final_heading/(norm(final_heading));
          headings(counter+1,:) = final_heading;
          next_movement_direction = final_heading;
end




%%%NOW PLOTTING PATH%%%
figure('Name','2D Path','NumberTitle','off');
scatter(positions(:,1),positions(:,2),[],[1:epochs+1]);
colormap(hsv);
colorbar;
hold on;
plot(positions(1,1),positions(1,2),'X','MarkerSize',25,'MarkerFaceColor', 'r');
line([-2 -2],[0 2.6],'color', 'k', 'LineStyle', '--');
line([-2 -1],[2.6 2.6],'color', 'k', 'LineStyle', '--');
line([-2 2],[1.6 1.6], 'color','k', 'LineStyle', '--');
line([-1 -1],[0 2.6], 'color','k', 'LineStyle', '--');
line([0 0],[0 1.6],'color', 'k', 'LineStyle', '--');
line([1 1],[0 1.6], 'color','k', 'LineStyle', '--');


%%%NOW PLOTTING HEADINGS%%%
figure('Name','2D Headings','NumberTitle','off');
quiver(positions(:,1),positions(:,2),headings(:,1),headings(:,2),0.4,'b');
hold on;
plot(positions(1,1),positions(1,2),'X','MarkerSize',25,'MarkerFaceColor', 'r');
line([-2 -2],[0 2.6],'color', 'k', 'LineStyle', '--');
line([-2 -1],[2.6 2.6],'color', 'k', 'LineStyle', '--');
line([-2 2],[1.6 1.6], 'color','k', 'LineStyle', '--');
line([-1 -1],[0 2.6], 'color','k', 'LineStyle', '--');
line([0 0],[0 1.6],'color', 'k', 'LineStyle', '--');
line([1 1],[0 1.6], 'color','k', 'LineStyle', '--');

%%%NOW TO RECONSTRUCT THE 3D PATH FROM THIS 2D DATA%%%



headings3d = zeros(epochs+1,3);
positions3d = zeros(epochs+1,3);



for idx = 1:epochs+1
    

    if(positions(idx,1)<-1)      %Transforming WEST wall
        if(positions(idx,2)<1.6) %Below top portion
            positions3d(idx,1) = -0.5;                  
            positions3d(idx,2) = -(positions(idx,1)+1.5);   
            positions3d(idx,3) = positions(idx,2); 
            
          
        else %top portion
            positions3d(idx,1) = positions(idx,2) - 2.1; 
            positions3d(idx,2) = -(positions(idx,1) + 1.5); 
            positions3d(idx,3) = 1.6;
         
            
        end
      
    elseif(positions(idx,1)<0) %Transforming SOUTH wall
            positions3d(idx,1) = positions(idx,1) + 0.5;
            positions3d(idx,2) = -0.5;
            positions3d(idx,3) = positions(idx,2);
            
         
            
            
    elseif(positions(idx,1)<1) %Transforming EAST wall
        positions3d(idx,1) = 0.5;
        positions3d(idx,2) = positions(idx,1) - 0.5;
        positions3d(idx,3) = positions(idx,2);
      
       
    else                       %Transforming NORTH wall
        positions3d(idx,1) = -(positions(idx,1) - 1.5);
        positions3d(idx,2) = 0.5;
        positions3d(idx,3) = positions(idx,2);
        
    end
    
        if(idx>1)
            headings3d(idx-1,:) = positions3d(idx,:) - positions3d(idx-1,:);
        end
end

headings3d(epochs+1,:) = headings3d(epochs,:); %bit artificial...

%Correcting headings to always lie in current plane - otherwise issues when
%moving between planes

for idx = 1:epochs+1
    if positions3d(idx,3) == 1.6
        headings3d(idx,3) = 0;
    end
    
    if positions3d(idx,2) == 0.5
        headings3d(idx,2) = 0;
    end
    
    if positions3d(idx,2) == -0.5
        headings3d(idx,2) = 0;
    end
    
     if positions3d(idx,1) == 0.5
        headings3d(idx,1) = 0;
    end
    
    if positions3d(idx,1) == -0.5
        headings3d(idx,1) = 0;
    end
end

%Normalising the headings too

for idx = 1:(epochs+1)
    headings3d(idx,:) = headings3d(idx,:)/(norm(headings3d(idx,:)));
end


%%%PLOTTING 3D PATH%%%
figure('Name','3D Path','NumberTitle','off')
%scatter3(positions3d(:,1),positions3d(:,2),positions3d(:,3),[],[1:epochs+1]);
%colormap(hsv);
%colorbar;
plot3(positions3d(:,1),positions3d(:,2),positions3d(:,3),'r');
xlabel('X axis');
ylabel('Y axis');

xlim([-1 1]);
ylim([-1 1]);
zlim([0 2]);

hold on
%Defining the cuboid corners
SouthWestBottom = [-0.5 -0.5   0];
SouthWestTop =    [-0.5 -0.5 1.6];
SouthEastBottom = [0.5  -0.5   0];
SouthEastTop =    [0.5  -0.5 1.6];


NorthWestBottom = [-0.5 0.5   0];
NorthWestTop =    [-0.5 0.5 1.6];
NorthEastBottom = [0.5  0.5   0];
NorthEastTop =    [0.5  0.5 1.6];

%Plotting south cube face
line([SouthWestBottom(1) SouthWestTop(1)],[SouthWestBottom(2) SouthWestTop(2)],[SouthWestBottom(3) SouthWestTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthEastBottom(1) SouthEastTop(1)],[SouthEastBottom(2) SouthEastTop(2)],[SouthEastBottom(3) SouthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthWestBottom(1) SouthEastBottom(1)],[SouthWestBottom(2) SouthEastBottom(2)],[SouthWestBottom(3) SouthEastBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthWestTop(1) SouthEastTop(1)],[SouthWestTop(2) SouthEastTop(2)],[SouthWestTop(3) SouthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');

%Plotting north cube face
line([NorthWestBottom(1) NorthWestTop(1)],[NorthWestBottom(2) NorthWestTop(2)],[NorthWestBottom(3) NorthWestTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([NorthEastBottom(1) NorthEastTop(1)],[NorthEastBottom(2) NorthEastTop(2)],[NorthEastBottom(3) NorthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([NorthWestBottom(1) NorthEastBottom(1)],[NorthWestBottom(2) NorthEastBottom(2)],[NorthWestBottom(3) NorthEastBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([NorthWestTop(1) NorthEastTop(1)],[NorthWestTop(2) NorthEastTop(2)],[NorthWestTop(3) NorthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');

%Plotting west cube face
line([SouthWestBottom(1) NorthWestBottom(1)],[SouthWestBottom(2) NorthWestBottom(2)],[SouthWestBottom(3) NorthWestBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthWestTop(1) NorthWestTop(1)],[SouthWestTop(2) NorthWestTop(2)],[SouthWestTop(3) NorthWestTop(3)]...
    ,'color', 'k', 'LineStyle', '--');

%Plotting east cube face
line([SouthEastBottom(1) NorthEastBottom(1)],[SouthEastBottom(2) NorthEastBottom(2)],[SouthEastBottom(3) NorthEastBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthEastTop(1) NorthEastTop(1)],[SouthEastTop(2) NorthEastTop(2)],[SouthEastTop(3) NorthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');


daspect([1 1 1]);
%%%PLOTTING 3D HEADINGS%%%
figure('Name','3D Headings','NumberTitle','off')
quiver3(positions3d(:,1),positions3d(:,2),positions3d(:,3),headings3d(:,1),headings3d(:,2),headings3d(:,3),0.4,'b');

xlabel('X axis');
ylabel('Y axis');

xlim([-1 1]);
ylim([-1 1]);
zlim([0 2]);

hold on

%Plotting south cube face
line([SouthWestBottom(1) SouthWestTop(1)],[SouthWestBottom(2) SouthWestTop(2)],[SouthWestBottom(3) SouthWestTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthEastBottom(1) SouthEastTop(1)],[SouthEastBottom(2) SouthEastTop(2)],[SouthEastBottom(3) SouthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthWestBottom(1) SouthEastBottom(1)],[SouthWestBottom(2) SouthEastBottom(2)],[SouthWestBottom(3) SouthEastBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthWestTop(1) SouthEastTop(1)],[SouthWestTop(2) SouthEastTop(2)],[SouthWestTop(3) SouthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');

%Plotting north cube face
line([NorthWestBottom(1) NorthWestTop(1)],[NorthWestBottom(2) NorthWestTop(2)],[NorthWestBottom(3) NorthWestTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([NorthEastBottom(1) NorthEastTop(1)],[NorthEastBottom(2) NorthEastTop(2)],[NorthEastBottom(3) NorthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([NorthWestBottom(1) NorthEastBottom(1)],[NorthWestBottom(2) NorthEastBottom(2)],[NorthWestBottom(3) NorthEastBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([NorthWestTop(1) NorthEastTop(1)],[NorthWestTop(2) NorthEastTop(2)],[NorthWestTop(3) NorthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');

%Plotting west cube face
line([SouthWestBottom(1) NorthWestBottom(1)],[SouthWestBottom(2) NorthWestBottom(2)],[SouthWestBottom(3) NorthWestBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthWestTop(1) NorthWestTop(1)],[SouthWestTop(2) NorthWestTop(2)],[SouthWestTop(3) NorthWestTop(3)]...
    ,'color', 'k', 'LineStyle', '--');

%Plotting east cube face
line([SouthEastBottom(1) NorthEastBottom(1)],[SouthEastBottom(2) NorthEastBottom(2)],[SouthEastBottom(3) NorthEastBottom(3)]...
    ,'color', 'k', 'LineStyle', '--');
line([SouthEastTop(1) NorthEastTop(1)],[SouthEastTop(2) NorthEastTop(2)],[SouthEastTop(3) NorthEastTop(3)]...
    ,'color', 'k', 'LineStyle', '--');


daspect([1 1 1]);
%Writing path data for reading into c code.

formatSpec = '%f\n';

fileID = fopen('pos2d_x.txt','w');
fprintf(fileID, formatSpec, positions(:,1));
fclose(fileID);

fileID = fopen('pos2d_y.txt','w');
fprintf(fileID, formatSpec, positions(:,2));
fclose(fileID);

fileID = fopen('positions_x_short.txt','w');
fprintf(fileID, formatSpec, positions3d(:,1));
fclose(fileID);

fileID = fopen('positions_y_short.txt','w');
fprintf(fileID, formatSpec, positions3d(:,2));
fclose(fileID);

fileID = fopen('positions_z_short.txt','w');
fprintf(fileID, formatSpec, positions3d(:,3));
fclose(fileID);

fileID = fopen('headings_x_short.txt','w');
fprintf(fileID, formatSpec, headings3d(:,1));
fclose(fileID);

fileID = fopen('headings_y_short.txt','w');
fprintf(fileID, formatSpec, headings3d(:,2));
fclose(fileID);

fileID = fopen('headings_z_short.txt','w');
fprintf(fileID, formatSpec, headings3d(:,3));
fclose(fileID);

%%Saving extended headings and path


extended_positions = zeros((epochs*100)+1,3);
extended_headings = zeros((epochs*100)+1,3);



for idx = 0:epochs-1
    for jdx = 1:100
        extended_positions((idx*100)+jdx,1) = positions3d(idx+1,1);
        extended_positions((idx*100)+jdx,2) = positions3d(idx+1,2);
        extended_positions((idx*100)+jdx,3) = positions3d(idx+1,3);
        
        extended_headings((idx*100)+jdx,1) = headings3d(idx+1,1);
        extended_headings((idx*100)+jdx,2) = headings3d(idx+1,2);
        extended_headings((idx*100)+jdx,3) = headings3d(idx+1,3);
    end
end

extended_positions((epochs*100)+1,1) = positions3d(epochs+1,1);
extended_positions((epochs*100)+1,2) = positions3d(epochs+1,2);
extended_positions((epochs*100)+1,3) = positions3d(epochs+1,3);

extended_headings((epochs*100)+1,1) = headings3d(epochs+1,1);
extended_headings((epochs*100)+1,2) = headings3d(epochs+1,2);
extended_headings((epochs*100)+1,3) = headings3d(epochs+1,3);


fileID = fopen('positions_x.txt','w');
fprintf(fileID, formatSpec, extended_positions(:,1));
fclose(fileID);

fileID = fopen('positions_y.txt','w');
fprintf(fileID, formatSpec, extended_positions(:,2));
fclose(fileID);

fileID = fopen('positions_z.txt','w');
fprintf(fileID, formatSpec, extended_positions(:,3));
fclose(fileID);

fileID = fopen('headings_x.txt','w');
fprintf(fileID, formatSpec, extended_headings(:,1));
fclose(fileID);

fileID = fopen('headings_y.txt','w');
fprintf(fileID, formatSpec, extended_headings(:,2));
fclose(fileID);

fileID = fopen('headings_z.txt','w');
fprintf(fileID, formatSpec, extended_headings(:,3));
fclose(fileID);

end
