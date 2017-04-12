function egocentric_yaw

head1 = [0,0.050000,0];
head2 = [0.015611,0.047359,-0.002368];

pos1 = [0,0,1];
pos2 = [0,0.049938,0.998752];

path_generated_rotation = 9.169444;

surface_normal1 = 2*pos1;

surface_normal2 = 2*pos2;

cross(surface_normal1, surface_normal2)
norm(cross(surface_normal1,surface_normal2))


        rot_angle = atan2(norm(cross(surface_normal1,surface_normal2)), dot(surface_normal1,surface_normal2));
        rot_axis = cross(surface_normal2,surface_normal1);
    
        rot_axis = rot_axis/norm(rot_axis);
    
        
        initial_quat = zeros(1,4);
        initial_quat(2) = head1(1);
        initial_quat(3) = head1(2);
        initial_quat(4) = head1(3);
        
        rotation_quat = zeros(1,4);
        rotation_quat(1) = cos(rot_angle/2);
        rotation_quat(2:4) = (sin(rot_angle/2) * rot_axis);
        
        %multiplying rotation quat by initial quat to get intermediate heading
        
        intermediate_quat = zeros(1,4);
        
        intermediate_quat(1) = (rotation_quat(1) * initial_quat(1)) - (rotation_quat(2) * initial_quat(2))...
            - (rotation_quat(3) * initial_quat(3)) - (rotation_quat(4) * initial_quat(4));
        
        intermediate_quat(2) = (rotation_quat(1) * initial_quat(2)) + (rotation_quat(2) * initial_quat(1))...
            - (rotation_quat(3) * initial_quat(4)) + (rotation_quat(4) * initial_quat(3));
        
        intermediate_quat(3) = (rotation_quat(1) * initial_quat(3)) + (rotation_quat(2) * initial_quat(4))...
            + (rotation_quat(3) * initial_quat(1)) - (rotation_quat(4) * initial_quat(2));
        
        intermediate_quat(4) = (rotation_quat(1) * initial_quat(4)) - (rotation_quat(2) * initial_quat(3))...
            + (rotation_quat(3) * initial_quat(2)) + (rotation_quat(4) * initial_quat(1));
        
        
        %finding inverse rotation_quat
        sum = 0;
        for index=1:4
            sum = sum + rotation_quat(index) * rotation_quat(index);
        end
        
        rotation_quat_magnitude = sqrt(sum);
        
        %rotation quat conjugate
        rotation_quat_conjugate = zeros(1,4);
        rotation_quat_conjugate(1) = rotation_quat(1);
        
        for index=2:4
            rotation_quat_conjugate(index) = -rotation_quat(index);
        end
        
        %Now using above to inverse rotation quaternion
        inverse_quat = zeros(1,4);
        
        for index=1:4
            inverse_quat(index) = rotation_quat_conjugate(index)/rotation_quat_magnitude;
        end
        
        %%multiplying intermediate_quat by inverse rotation quat to get final_quat
        
        final_quat = zeros(1,4);
        
        
        final_quat(1) = (intermediate_quat(1) * inverse_quat(1)) - (intermediate_quat(2) * inverse_quat(2))...
            - (intermediate_quat(3) * inverse_quat(3)) - (intermediate_quat(4) * inverse_quat(4));
        
        final_quat(2) = (intermediate_quat(1) * inverse_quat(2)) + (intermediate_quat(2) * inverse_quat(1))...
            - (intermediate_quat(3) * inverse_quat(4)) + (intermediate_quat(4) * inverse_quat(3));
        
        final_quat(3) = (intermediate_quat(1) * inverse_quat(3)) + (intermediate_quat(2) * inverse_quat(4))...
            + (intermediate_quat(3) * inverse_quat(1)) - (intermediate_quat(4) * inverse_quat(2));
        
        final_quat(4) = (intermediate_quat(1) * inverse_quat(4)) - (intermediate_quat(2) * inverse_quat(3))...
            + (intermediate_quat(3) * inverse_quat(2)) + (intermediate_quat(4) * inverse_quat(1));
        
        
        rotated_reference = [final_quat(2),final_quat(3),final_quat(4)];
        
        
        figure()
        quiver3(0,0,0,head1(1),head1(2),head1(3),'b');

        hold on
        quiver3(0,0,0,head2(1),head2(2),head2(3),'b');
        
        quiver3(0,0,0,rotated_reference(1),rotated_reference(2),rotated_reference(3),'g');

        
        egocentric_rotation = atan2(norm(cross(rotated_reference,head2)), dot(rotated_reference,head2)) *(180/pi)
        
        egocentric_axis = cross(head2,rotated_reference);
        
        quiver3(0,0,0,egocentric_axis(1),egocentric_axis(2),egocentric_axis(3),'k');
        
        quiver3(0,0,0,surface_normal2(1)*0.01,surface_normal2(2)*0.01,surface_normal2(3)*0.01,'k');
        
        quiver3(0,0,0,surface_normal1(1)*0.01,surface_normal1(2)*0.01,surface_normal1(3)*0.01,'k');
        
        
end