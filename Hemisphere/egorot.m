idx = 4140

figure()
quiver3(0,0,0,0,0,1,'g')
hold on
xlabel('x axis');
ylabel('y axis');
quiver3(pos_x(idx),pos_y(idx),pos_z(idx),refvecs(idx,1),refvecs(idx,2),refvecs(idx,3),'b')
quiver3(pos_x(idx),pos_y(idx),pos_z(idx),head_x(idx),head_y(idx),head_z(idx),'k')
quiver3(pos_x(idx),pos_y(idx),pos_z(idx),pos_x(idx),pos_y(idx),pos_z(idx),'g');


quiver3(pos_x(idx+1),pos_y(idx+1),pos_z(idx+1),refvecs(idx+1,1),refvecs(idx+1,2),refvecs(idx+1,3),'r')
quiver3(pos_x(idx+1),pos_y(idx+1),pos_z(idx+1),head_x(idx+1),head_y(idx+1),head_z(idx+1),'k')
quiver3(pos_x(idx+1),pos_y(idx+1),pos_z(idx+1),pos_x(idx+1),pos_y(idx+1),pos_z(idx+1),'g');

refvec1 = [refvecs(idx,1) refvecs(idx,2) refvecs(idx,3)];
refvec1 = refvec1/(norm(refvec1));

dotidx = dot(refvec1,[pos_x(idx) pos_y(idx) pos_z(idx)]);

refvec2 = [refvecs(idx+1,1) refvecs(idx+1,2) refvecs(idx+1,3)];
refvec2 = refvec2/(norm(refvec2));

dotjdx = dot(refvec2,[pos_x(idx+1) pos_y(idx+1) pos_z(idx+1)]);

head1 = [head_x(idx) head_y(idx) head_z(idx)];
head2 = [head_x(idx+1) head_y(idx+1) head_z(idx+1)];

pos1 = [pos_x(idx) pos_y(idx) pos_z(idx)];
pos2 = [pos_x(idx+1) pos_y(idx+1) pos_z(idx+1)];

pos1 = pos1/(norm(pos1));
pos2 = pos2/(norm(pos2));

figure()

quiver3(0,0,0,head_x(idx),head_y(idx),head_z(idx),'b')
hold on
quiver3(0,0,0,refvecs(idx,1),refvecs(idx,2),refvecs(idx,3),'b')
quiver3(0,0,0,pos1(1),pos1(2),pos1(3),'color','b','LineStyle','--')

quiver3(0,0,0,head_x(idx+1),head_y(idx+1),head_z(idx+1),'r')
quiver3(0,0,0,refvecs(idx+1,1),refvecs(idx+1,2),refvecs(idx+1,3),'r')
quiver3(0,0,0,pos2(1),pos2(2),pos2(3),'color','r','LineStyle','--')
